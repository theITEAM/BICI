// Implements an MCMC algorithm

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
 
using namespace std;

#include "mcmc.hh"
#include "state.hh"
#include "utils.hh"

MCMC::MCMC(const Model &model, Output &output) : model(model), output(output), state(model)
{	
	const auto &de = model.details;
	nsample = de.sample;
	nburnin = (unsigned int)(0.2*nsample);
	thin_param = de.thinparam;
	thin_state = de.thinstate;
}
 
 
/// Runs the inference algorithm 
void MCMC::run()
{
	param_prop_init();
	
	//output.print_param(param_val);
	
	Particle part;
	
	auto like_obs_max = -LARGE;
	
	auto loop_max = 50u; 
	for(auto loop = 0u; loop < loop_max; loop++){ 
		auto initc_val = model.initc_sample();
		auto param_val = model.param_sample();
	
		state.simulate(param_val,initc_val);
	
		state.resample_ind();   // Resamples individual such that fixed events become correct
		
		if(state.like.obs > like_obs_max){
			like_obs_max = state.like.obs;
			part = state.generate_particle();
		}
	}
	
	state.set_particle(part);
	
	state.check("begin");
	
	//state.species[0].generate_A();
	
	long time_start = -clock();

	for(auto s = 0u; s < nsample; s++){
		if(s%10 == 0 && com_op) progress(s,nsample);		
		else{
			if(s%thin_param == 0) cout << s << " " << nsample << " sample\n";
		}
		
		if(s < nburnin){ burnin = true; if(s < nburnin/2) burnfac = 10; else burnfac = 1;}
		else burnin = false;
		
		if(s%thin_param == 0) output.param_sample(s,state);
		
		if(s%thin_state == 0) output.state_sample(s,state);
		
		update();
		
		if(burnin){
			for(auto &pro : proposal){
				if(pro.type == PARAM_PROP){
					pro.add_sample(state.param_val);
					if(s%10 == 0) pro.update_sampler();
				}
			}
		}
		
		if(s%10 == 0) state.check("end");
	}
	
	//state.species[0].check_transnum();
	
	diagnostics(clock()-time_start);
	
	//output.print_param(param_val);
}


/// Initialises some proposal distributions
void MCMC::param_prop_init()
{
	for(auto i = 0u; i < model.nparam_vec; i++){ // Univariate distributions
		const auto &pv = model.param_vec[i];
		const auto &par = model.param[pv.th];
		if(par.variety == PRIOR_PARAM || par.variety == DIST_PARAM){
			if(par.prior[pv.index].type != FIX_PR){
				vector <unsigned int> vec;
				vec.push_back(i);
	
				Proposal pp(PARAM_PROP,vec,model,output,1,burnfac);
				proposal.push_back(pp);
				
				/// Looks at adding MBPs
				vector <unsigned int> list;
				for(const auto &al : pp.affect_like){
					if(al.type == MARKOV_POP_AFFECT){
						add_to_vec(list,al.num);
					}
				}
			
				for(auto k = 0u; k < list.size(); k++){ // Adds one MBP per species affected
					auto pp2 = pp;
					pp2.p_prop = list[k];
					
					auto j = 0u;                          // Removes MARKOV_POP_AFFECT from list
					while(j < pp2.affect_like.size()){
						auto al = pp2.affect_like[j];
						if(al.type == MARKOV_POP_AFFECT && al.num == list[k]){
							pp2.affect_like.erase(pp2.affect_like.begin()+j);
						}
						else j++;
					}
		
					AffectLike al; al.type = LIKE_OBS_AFFECT; al.num = pp2.p_prop; al.num2 = UNSET;
					param_vec_add_affect(pp2.affect_like,al);
					
					pp2.type = MBP_PROP;
					pp2.name = "MBP "+pp2.name;
					pp2.mbp_population_affect();
					proposal.push_back(pp2);
				}
			}
		}
	}
	
	/// Does proposals which simultaneously change variances and ies
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		
		for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
			auto &ieg = sp.ind_eff_group[i];
			for(auto j = 0u; j < ieg.list.size(); j++){
				auto &par = model.param[ieg.omega[j][j]];
				if(par.variety != CONST_PARAM){
					vector <unsigned int> vec;
					vec.push_back(par.param_vec_ref[0]);
					vec.push_back(p);
					vec.push_back(ieg.list[j].index );
					
					Proposal pp(IE_VAR_PROP,vec,model,output,1,burnfac);
					proposal.push_back(pp);
				}
			}
		}
	}
	
	for(auto i = 0u; i < 10; i++){
		auto param_val = model.param_sample();
		for(auto &pro : proposal) pro.add_sample(param_val);
	}
	
	for(auto &pro : proposal) pro.update_sampler();
	
	
	for(auto loop = 0u; loop < 1; loop++){
		for(auto p = 0u; p < model.nspecies; p++){
			if(model.species[p].type == POPULATION){
				vector <unsigned int> vec;
				vec.push_back(p);
				Proposal pp(MBPII_PROP,vec,model,output,1,burnfac);
				
				AffectLike al; al.type = LIKE_OBS_AFFECT; al.num = p; al.num2 = UNSET;
				param_vec_add_affect(pp.affect_like,al);
				
				proposal.push_back(pp);
			}
		}
	}
				
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			vector <unsigned int> vec; vec.push_back(p);
			Proposal pp(IND_PROP,vec,model,output,0.3,burnfac);
			proposal.push_back(pp);
			
			Proposal pp2(INDSAMP_PROP,vec,model,output,0.1,burnfac);
			proposal.push_back(pp2);
		}
	}
	
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.ind_effect.size(); i++){
				vector <unsigned int> vec;
				vec.push_back(p); vec.push_back(i);
				Proposal pp(IE_PROP,vec,model,output,0.3,burnfac);
				proposal.push_back(pp);
			}
		}
	}
	
	for(auto &pro : proposal) pro.print_info();

	//emsg("done"); //zz
}


/// Updates all the proposals
void MCMC::update()
{
	for(auto &pro : proposal){
		//cout << pro.name << " proposal\n";
		pro.update(state,burnin);
	}
}



/// Outputs diagnostics about the proposals
void MCMC::diagnostics(long total_time) const
{                        
	auto file = output.outputdir+"/diagnostics.txt";
	ofstream fout(file);
	output.check_open(fout,file);
	
	fout << "MCMC PROPOSAL DIAGNOSTICS" << endl << endl;
	
	for(auto &pro : proposal){
		fout << pro.diagnostics(total_time) << endl;
	}
	
	fout <<  "Check CPU  time: " << cpu_percent(state.check_timer,total_time) << endl;
	fout << endl;

	fout <<  "Total CPU  time: " << total_time/CLOCKS_PER_SEC << " seconds" << endl;
}
