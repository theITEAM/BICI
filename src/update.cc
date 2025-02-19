// Implements an update to the state making use of MCMC proposals

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm> 
 
using namespace std;

#include "update.hh"
#include "state.hh"
#include "utils.hh"
#include "matrix.hh"

Update::Update(const Model &model, const BurnInfo &burn_info, Output &output) : model(model), burn_info(burn_info), output(output)
{	
}


/// Updates all the proposals
void Update::run(unsigned int s, State &state)
{
	state.sample = s;
 
	auto pl = false;
	
	for(auto &pro : proposal){ 
		if(pro.on){
			if(pl) cout << pro.name << " " <<  state.sample << " proposal" << endl;
			
			if(true){
				pro.update(state);
			}
			else{
				switch(pro.type){
				case PAR_EVENT_FORWARD_PROP:
				case PAR_EVENT_FORWARD_SQ_PROP:
				case PAR_EVENT_BACKWARD_SQ_PROP:
				//case IND_ADD_REM_PROP:
					pro.update(state);
					break;
			
				default:
					break;
				}
			}
		
			if(pl) state.check(" After prop check");
		}
	}
	
	if(testing == true && s%CHECK_THIN == CHECK_THIN-1) state.check("end");
}


/// Initialises some proposal distributions
void Update::init()
{	
	cor_matrix.init(model.nparam_vec);

	for(auto i = 0u; i < model.nparam_vec; i++){     // Univariate distributions
		if(model.param_vec[i].prop_pos){
			vector <unsigned int> vec;
			vec.push_back(i);
			add_parameter_prop(vec);
		}
	}
	
	// Proposals on individual effects
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			for(auto i = 0u; i < sp.ind_effect.size(); i++){
				vector <unsigned int> vec;
				vec.push_back(p); vec.push_back(i);
				Proposal pp(IE_PROP,vec,model,output,0.3,burn_info);
				proposal.push_back(pp);
			}
		}
	}
	
	// Does proposals which simultaneously change variances and ies
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
					vec.push_back(ieg.list[j].index);
					
					Proposal pp(IE_VAR_PROP,vec,model,output,1,burn_info);
					proposal.push_back(pp);
				}
			}
			
			if(ieg.list.size() > 0){
				for(auto j = 0u; j < ieg.list.size(); j++){
					for(auto k = j+1; k < ieg.list.size(); k++){
						auto &par = model.param[ieg.omega[j][k]];
						if(par.variety != CONST_PARAM){
							vector <unsigned int> vec;
							vec.push_back(par.param_vec_ref[0]);
							vec.push_back(p);
							vec.push_back(i);
							vec.push_back(j);
							vec.push_back(k);		
							Proposal pp(IE_COVAR_PROP,vec,model,output,1,burn_info);
							proposal.push_back(pp);
						}
					}
				}
			}
		}
	}
	
	// Adds MBPs which don't change parameters
	for(auto loop = 0u; loop < 1; loop++){
		for(auto p = 0u; p < model.nspecies; p++){
			if(model.species[p].type == POPULATION){
				vector <unsigned int> vec;
				vec.push_back(p);
				Proposal pp(MBPII_PROP,vec,model,output,PROP_MBPII_W,burn_info);
				model.add_like_obs_affect(p,pp.affect_like);
				proposal.push_back(pp);
			}
		}
	}
	
	// Adds MBPs which change initial conditions
	for(auto loop = 0u; loop < 1; loop++){
		for(auto p = 0u; p < model.nspecies; p++){
			const auto &sp = model.species[p];
			const auto &ic = sp.init_cond;
			if(sp.type == POPULATION && ic.type == INIT_POP_DIST){
				if(sp.ncla > 1){
					vector <unsigned int> vec;
					vec.push_back(p);
					Proposal pp(MBP_IC_RESAMP_PROP,vec,model,output,PROP_MBPII_W,burn_info);		
					model.add_like_obs_affect(p,pp.affect_like);
					proposal.push_back(pp);
				}					
					
				auto foc_cl = ic.focal_cl;
				if(foc_cl == UNSET){
					vector <unsigned int> vec;
					vec.push_back(p);
					Proposal pp(MBP_IC_POPTOTAL_PROP,vec,model,output,PROP_MBPIC_W,burn_info);
					AffectLike al; al.type = PRIOR_INIT_COND_AFFECT; al.num = UNSET; al.num2 = UNSET;
					param_vec_add_affect(pp.affect_like,al);		
					model.add_like_obs_affect(p,pp.affect_like);
					proposal.push_back(pp);
				}
				else{
					const auto &claa = sp.cla[foc_cl];
					for(auto c = 0u; c < claa.comp.size(); c++){
						if(!claa.comp[c].erlang_hidden){
							vector <unsigned int> vec;
							vec.push_back(p); vec.push_back(foc_cl); vec.push_back(c);
							
							Proposal pp(MBP_IC_POP_PROP,vec,model,output,PROP_MBPIC_W,burn_info);
							AffectLike al; al.type = PRIOR_INIT_COND_AFFECT; al.num = UNSET; al.num2 = UNSET;
							param_vec_add_affect(pp.affect_like,al);		
							model.add_like_obs_affect(p,pp.affect_like);
							proposal.push_back(pp);
						}
					}
				}
			}
		}
	}
	
	// Makes proposals which make local changes to population-based model
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		if(sp.type == POPULATION){
			// Finds potential local transition combinations
			auto popN_fixed = sp.is_pop_num_fixed();
		
			if(!popN_fixed){
				vector <unsigned int> vec; vec.push_back(p);
				Proposal pp(POP_IC_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
				proposal.push_back(pp);
			}
			
			if(sp.tra_gl.size() > 0){
				{
					vector <unsigned int> vec; vec.push_back(p); 
					Proposal pp(POP_SINGLE_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp);
				}
				
				{
					auto list = sp.find_add_rem_list();
					for(auto &li : list){
						auto vec = li.tr_list;
						vec.push_back(p);
						Proposal pp(POP_ADD_REM_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
						proposal.push_back(pp);
					}
				}

				if(sp.init_cond.type == INIT_POP_DIST){	
					vector <unsigned int> vec; vec.push_back(p);
					Proposal pp(POP_IC_SWAP_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp);
				}
				
				{
					vector <unsigned int> vec; vec.push_back(p);
					Proposal pp(POP_MOVE_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp);
				}
				
				{
					vector <unsigned int> vec; vec.push_back(p);
					Proposal pp(POP_END_LOCAL_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
					proposal.push_back(pp);
				}
			}
		}
	}
	
	// Makes poposals in init cont which change initial condition fractions
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		const auto &ic = sp.init_cond;
		if(ic.type == INIT_POP_DIST){
			vector <unsigned int> vec; vec.push_back(p);
			Proposal pp(INIT_COND_FRAC_PROP,vec,model,output,PROP_LOCAL_W,burn_info);
			proposal.push_back(pp);
		}
	}
	
	if(true){
		for(auto p = 0u; p < model.nspecies; p++){
			auto &sp = model.species[p];
			auto st = model.samp_type;
			
			if(sp.type == INDIVIDUAL){	
				vector <unsigned int> vec; vec.push_back(p);
				{				
					Proposal pp(IND_EVENT_TIME_PROP,vec,model,output,1,burn_info);
					//if(st == LOCAL_SAMP || st == ALL_SAMP) proposal.push_back(pp); 
					proposal.push_back(pp); 
				}
				
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(sp.cla[cl].swap_rep.size() > 0){
						vector <unsigned int> vec2; vec2.push_back(p); vec2.push_back(cl);
						Proposal pp(IND_LOCAL_PROP,vec2,model,output,1,burn_info);
						if(st == LOCAL_SAMP || st == ALL_SAMP) proposal.push_back(pp); 
					}
				}	
				
				if(sp.nindividual_in > 0){	 
					Proposal pp2(IND_OBS_SAMP_PROP,vec,model,output,0.1,burn_info);
					if(st == SAMP_SAMP || st == ALL_SAMP) proposal.push_back(pp2);
					
					Proposal pp3(IND_OBS_RESIM_PROP,vec,model,output,0.1,burn_info);
					if(st == SIM_SAMP || st == ALL_SAMP) proposal.push_back(pp3);
					
					Proposal pp4(IND_OBS_RESIM_SINGLE_PROP,vec,model,output,0.1,burn_info);
					if(st == SIM_CL_SAMP || st == ALL_SAMP) proposal.push_back(pp4);
				}
				
				{
					Proposal pp(IND_UNOBS_RESIM_PROP,vec,model,output,1,burn_info);
					if(st == SIM_SAMP || st == ALL_SAMP) proposal.push_back(pp);
				}
				
				if(sp.contains_source || sp.init_cond.type == INIT_POP_DIST){
					if(sp.trans_tree){ // For trans-tree only add/rem a single individual
						Proposal pp(IND_ADD_REM_TT_PROP,vec,model,output,1,burn_info);
						proposal.push_back(pp);
					}
					else{              // Otherwise add/rem multiple individuals
						Proposal pp(IND_ADD_REM_PROP,vec,model,output,1,burn_info);
						proposal.push_back(pp);
					} 
				}
				
				if(sp.contains_source && sp.init_cond.type != INIT_POP_FIXED){
					Proposal pp(IND_OBS_SWITCH_ENTER_SOURCE_PROP,vec,model,output,1,burn_info);
					proposal.push_back(pp);
				}
			}
		}
	}

	if(false){
		// Joint proposals on events and parameters
		for(auto p = 0u; p < model.nspecies; p++){
			auto &sp = model.species[p];
			if(sp.type == INDIVIDUAL){
				for(const auto &pej : sp.par_event_joint){
					vector <unsigned int> vec;
					vec.push_back(p);
					vec.push_back(pej.th);
					for(auto tr : pej.tr_list) vec.push_back(tr);
					switch(pej.dir){
					case FORWARD:
						{
							Proposal pp(PAR_EVENT_FORWARD_PROP,vec,model,output,0.3,burn_info);
							//proposal.push_back(pp);
						}
						break;
						
					case FORWARD_SQ:
						{
							Proposal pp(PAR_EVENT_FORWARD_SQ_PROP,vec,model,output,0.3,burn_info);
							//proposal.push_back(pp);
						}
						break;
						
					case BACKWARD_SQ:
						{
							Proposal pp(PAR_EVENT_BACKWARD_SQ_PROP,vec,model,output,0.3,burn_info);
							//proposal.push_back(pp);
						}
						break;
					}
				}
			}
		}		
	}
	
	// Adds a transmission tree proposal
	if(model.trans_tree){
		vector <unsigned int> vec;
		{
			Proposal pp(TRANS_TREE_PROP,vec,model,output,1,burn_info);
			proposal.push_back(pp);
		}
		
		{
			Proposal pp(TRANS_TREE_SWAP_INF_PROP,vec,model,output,1,burn_info);
			proposal.push_back(pp);
		}
		
		if(model.genetic_data.on){
			{
				Proposal pp(TRANS_TREE_MUT_PROP,vec,model,output,1,burn_info);
				proposal.push_back(pp);
			}
			
			{
				Proposal pp(TRANS_TREE_MUT_LOCAL_PROP,vec,model,output,1,burn_info);
				proposal.push_back(pp);
			}
		}
	}
		
	// Adds proposals to correct observed transitions
	for(auto p = 0u; p < model.nspecies; p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			auto fl = false;
			for(auto i = 0u; i < sp.nindividual_obs; i++){
				for(const auto &ob : sp.individual[i].obs){
					switch(ob.type){
					case OBS_SINK_EV: case OBS_TRANS_EV: fl = true; break;
					default: break;
					}
				}
			}
			
			if(fl){			
				vector <unsigned int> vec; vec.push_back(p);
				Proposal pp(CORRECT_OBS_TRANS_PROP,vec,model,output,UNSET,burn_info);
				proposal.push_back(pp);
			}
		}
	}				
					
	for(auto i = 0u; i < 10; i++){
		auto param_val = model.param_sample();
		cor_matrix.add_sample(param_val);
	}
		
	for(auto &pro : proposal) pro.update_sampler(cor_matrix);
	
	for(auto &pro : proposal) pro.calculate_affect_spline(); 
	
	if(1 == 0) cout << "print_info_turned off in code" << endl;
	else{
		string te = "";
		for(auto &pro : proposal) te += pro.print_info();
		output.prop_summary(te);
	}
}


/// Adds a parameter 
void Update::add_parameter_prop(const vector <unsigned int> &vec)
{
	Proposal pp(PARAM_PROP,vec,model,output,1,burn_info);
	proposal.push_back(pp);
	
	/// Looks at adding MBPs
	vector <unsigned int> list;
	for(const auto &al : pp.affect_like){
		if(al.type == MARKOV_POP_AFFECT){
			add_to_vec(list,al.num);
		}
	}

	for(auto k = 0u; k < list.size(); k++){    // Adds one MBP per species affected
		auto pp2 = pp;
		pp2.p_prop = list[k];
		
		auto j = 0u;                             // Removes MARKOV_POP_AFFECT from list
		while(j < pp2.affect_like.size()){
			auto al = pp2.affect_like[j];
			if(al.type == MARKOV_POP_AFFECT && al.num == list[k]){
				pp2.affect_like.erase(pp2.affect_like.begin()+j);
			}
			else j++;
		}

		model.add_like_obs_affect(pp2.p_prop,pp2.affect_like);
	
		pp2.type = MBP_PROP;
		pp2.name = "MBP "+pp2.name;
		pp2.mbp_population_affect();
		proposal.push_back(pp2);
	}
}


/// Outputs diagnostics about the proposals
void Update::diagnostics(long total_time, const State &state) const
{                        
	auto file = output.outputdir+"/diagnostics.txt";
	ofstream fout(file);
	output.check_open(fout,file);

	{	
		fout << "MCMC PROPOSAL DIAGNOSTICS" << endl << endl;

		auto sum_samp = 0.0, sum_prop = 0.0, sum_update = 0.0, sum_mh = 0.0;
		auto sum_trans_tree = 0.0, sum_trans_tree_swap_inf = 0.0;
		auto sum_trans_tree_mut = 0.0, sum_trans_tree_mut_local = 0.0;
		auto sum_init_ind_local = 0.0, sum_ind_local = 0.0;
		auto sum_temp1 = 0.0, sum_temp2 = 0.0, sum_temp3 = 0.0, sum_temp4 = 0.0, sum_temp5 = 0.0;
		auto p_sel = UNSET;
		for(auto &pro : proposal){
			sum_prop += pro.timer[PROP_TIMER];
			sum_samp += pro.timer[SAMPLE_TIMER];
			sum_update += pro.timer[UPDATESAMP_TIMER];
			sum_mh += pro.timer[MH_TIMER];
			sum_init_ind_local += pro.timer[INIT_IND_LOCAL_TIMER];
			sum_ind_local += pro.timer[IND_LOCAL_TIMER];
			sum_trans_tree += pro.timer[TRANS_TREE_TIMER];
			sum_trans_tree_swap_inf += pro.timer[TRANS_TREE_SWAP_INF_TIMER];
			sum_trans_tree_mut += pro.timer[TRANS_TREE_MUT_TIMER];
			sum_trans_tree_mut_local += pro.timer[TRANS_TREE_MUT_LOCAL_TIMER];
			sum_temp1 += pro.timer[TEMP1_TIMER];
			sum_temp2 += pro.timer[TEMP2_TIMER];
			sum_temp3 += pro.timer[TEMP3_TIMER];
			sum_temp4 += pro.timer[TEMP4_TIMER];
			sum_temp5 += pro.timer[TEMP5_TIMER];
			
			auto p = pro.p_prop;
			if(p != UNSET && p != p_sel && model.species.size() > 1){
				fout << "<< For species: " << model.species[p].name << " >>" << endl;
				p_sel = p;
			}
			
			fout << pro.diagnostics(total_time);
			
			fout << endl;
		}
		fout << "Overall: proposal - " << cpu_percent(sum_prop,total_time) << "  ";
		fout << "sample - " << cpu_percent(sum_samp,total_time) 
				 << "  update - " << cpu_percent(sum_update,total_time)
				 << "  mh - " << cpu_percent(sum_mh,total_time)
				 << "  init ind local - " << cpu_percent(sum_init_ind_local,total_time)
				 << "  ind local after filt - " << cpu_percent(sum_ind_local,total_time)
				 << "  trans-tree - " << cpu_percent(sum_trans_tree,total_time)
				 << "  trans-tree-mut - " << cpu_percent(sum_trans_tree_mut,total_time)
				 << "  trans-tree-mut-local - " << cpu_percent(sum_trans_tree_mut_local,total_time)
				 << "  temp1 - " << cpu_percent(sum_temp1,total_time)
				 << "  temp2 - " << cpu_percent(sum_temp2,total_time)
				 << "  temp3 - " << cpu_percent(sum_temp3,total_time)
				 << "  temp4 - " << cpu_percent(sum_temp4,total_time)
				 << "  temp5 - " << cpu_percent(sum_temp5,total_time)
				 << endl;
		fout << endl;	
	}
	
	{
		fout << "UPDATE AND RESTORE DIAGNOSTICS" << endl << endl;
		auto sum_update = 0.0, sum_restore = 0.0;
		for(auto i = 0u; i < AFFECT_MAX; i++){	
			switch(i){
			case SPLINE_PRIOR_AFFECT: fout << "Spline prior"; break; 
			case PRIOR_AFFECT: fout << "Prior"; break; 
			case DIST_AFFECT: fout << "Distribution"; break; 
			case SPLINE_AFFECT: fout << "Spline"; break;  
			case EXP_FE_AFFECT: fout << "Exp FE"; break; 
			case DIV_VALUE_AFFECT: fout << "Div value"; break; 
			case DIV_VALUE_FAST_AFFECT: fout << "Div value fast"; break; 
			case DIV_VALUE_LINEAR_AFFECT: fout << "Div value linear"; break; 
			case MARKOV_LIKE_AFFECT: fout << "Markov like"; break; 
			case POP_AFFECT: fout << "Population"; break; 
			case NM_TRANS_AFFECT: fout << "NM trans"; break; 
			case NM_TRANS_BP_AFFECT: fout << "NM trans bp"; break; 
			case NM_TRANS_INCOMP_AFFECT: fout << "NM trans incomp"; break; 
			case OMEGA_AFFECT: fout << "Omega"; break; 
			case EXP_IE_AFFECT: fout << "Exp IE"; break; 
			case LIKE_IE_AFFECT: fout << "Like IE"; break; 
			case INDFAC_INT_AFFECT: fout << "Individual factor"; break; 
			case MARKOV_POP_AFFECT: fout << "Markov Pop"; break; 
			case OBS_EQN_AFFECT: fout << "Obs Eqn"; break; 
			case LIKE_OBS_IND_AFFECT: fout << "Like Obs Ind"; break; 
			case LIKE_OBS_POP_AFFECT: fout << "Like Obs Pop"; break; 
			case LIKE_OBS_POP_TRANS_AFFECT: fout << "Like Obs Pop Trans"; break; 
			case LIKE_UNOBS_TRANS_AFFECT: fout << "Like unobs trans"; break;
			case POP_DATA_CGL_TGL_AFFECT: fout << "Pop data cgl tgl"; break;
			case LIKE_INIT_COND_AFFECT: fout << "Like init cond"; break;
			case PRIOR_INIT_COND_AFFECT: fout << "prior init cond"; break;
			case LIKE_GENETIC_PROCESS_AFFECT: fout << "Like genetic process"; break; 
			case GENETIC_VALUE_AFFECT: fout << "Genetic value"; break; 
			case LIKE_GENETIC_OBS_AFFECT: fout << "Like genetic obs"; break; 
			}
			fout << "  ";
			fout << "update - " << cpu_percent(state.update_timer[i],total_time) << "  ";
	
			fout << "restore - " << cpu_percent(state.restore_timer[i],total_time) << endl;
			
			sum_update += state.update_timer[i];
			sum_restore += state.restore_timer[i];
		}
		fout << "Overall: update - " << cpu_percent(sum_update ,total_time) << "  ";
		fout << "restore - " << cpu_percent(sum_restore,total_time) << endl;
			
		fout << endl;
	}

	for(const auto &ssp : state.species){
		fout << "Linear " << cpu_percent(ssp.timer[LINEAR_TIMER],total_time) << endl;
		fout << "Linear zero " << cpu_percent(ssp.timer[LINEARZERO_TIMER],total_time) << endl;
		fout << "Splinechange " << cpu_percent(ssp.timer[SPLINECHANGE_TIMER],total_time) << endl;
		fout << "Splinechangegrad " << cpu_percent(ssp.timer[SPLINECHANGEGRAD_TIMER],total_time) << endl;
	}
	
	fout << "Output param " << cpu_percent(output.timer[PARAM_OUTPUT],total_time) << endl;
	fout << "Output state " << cpu_percent(output.timer[STATE_OUTPUT],total_time) << endl;
	
	fout <<  "Gen dif CPU  time: " << cpu_percent(state.timer[GEN_DIF_TIMER],total_time) << endl;
	fout <<  "trans tree like CPU  time: " << cpu_percent(state.timer[TRANS_TREE_LIKE_TIMER],total_time) << endl;
	fout <<  "genetic proc like CPU  time: " << cpu_percent(state.timer[GENETIC_PROC_LIKE_TIMER],total_time) << endl;
	fout <<  "genetic obs like CPU  time: " << cpu_percent(state.timer[GENETIC_OBS_LIKE_TIMER],total_time) << endl;
	
	fout <<  "ind update: " << cpu_percent(state.timer[IND_TIMER],total_time) << endl;
	
	fout <<  "ind pop update: " << cpu_percent(state.timer[IND_POP_UPDATE_TIMER],total_time) << endl;
	
	fout <<  "Update sampler: " << cpu_percent(state.timer[UPDATE_SAMPLER_TIMER],total_time) << endl;

	fout <<  "Restore ind: " << cpu_percent(state.timer[RESTORE_TIMER],total_time) << endl;
	fout << endl;
	
	fout <<  "Check CPU  time: " << cpu_percent(state.timer[CHECK_TIMER],total_time) << endl;
	fout << endl;

	fout <<  "Total CPU  time: " << total_time/CLOCKS_PER_SEC << " seconds" << endl;
	
	if(false){
		for(const auto &ssp : state.species){
			fout << ssp.source_sampler.print();
		}
	}
}


/// Used to order proposal speeds
bool PS_ord (PropSpeed ps1, PropSpeed ps2)                      
{ return (ps1.time_per_prop < ps2.time_per_prop); };  


/// Sets the proposal probabilities based on the amount of CPU time taken 
void Update::set_proposal_prob()
{
	vector <PropSpeed> list;
	for(auto i = 0u; i < proposal.size(); i++){
		const auto &pro = proposal[i];
		if(pro.number != 0){
			PropSpeed ps; ps.i = i; ps.time_per_prop = pro.timer[PROP_TIMER]/pro.number;
			list.push_back(ps);
		}
	}
	
	sort(list.begin(),list.end(),PS_ord);

	// List all possible transitions which could have reduced probabilities
	vector <PropType> pos = {POP_ADD_REM_LOCAL_PROP,POP_MOVE_LOCAL_PROP,POP_IC_LOCAL_PROP, POP_END_LOCAL_PROP, POP_SINGLE_LOCAL_PROP, POP_IC_PROP, POP_IC_SWAP_PROP,IND_OBS_SAMP_PROP,IND_OBS_RESIM_PROP};
	
	auto jmin = list.size();

	while(jmin > 0){
		if(list[jmin-1].time_per_prop == 0) break;
		
		auto ty = proposal[list[jmin-1].i].type;
		auto kmax = pos.size();
		auto k = 0u; while(k < kmax && pos[k] != ty) k++;
		if(k == kmax) break;
		jmin--;
	}

	for(auto j = 0u; j < jmin; j++){
		auto &pro = proposal[list[j].i];
		pro.prop_prob = pro.prop_weight;
	}
	
	if(jmin < list.size()){
		auto ref = list[jmin].time_per_prop;
		for(auto j = jmin; j < list.size(); j++){
			auto &pro = proposal[list[j].i];
			auto prob = (ref/list[j].time_per_prop)*pro.prop_weight;
			if(prob < PROPOSAL_PROB_LIM) prob = PROPOSAL_PROB_LIM;
			pro.prop_prob = prob;
		}
	}
}


/// Checks if new proposals are added which join parameters
void Update::check_join_proposal()
{
	auto M = cor_matrix.calculate_cor_matrix();

	print_matrix("mat",M);

	auto f = exp(-(cor_matrix.get_n()/400.0));
	auto thresh = f + (1-f)*PROP_JOIN_COR_MIN;

	auto N = M.size();
	for(auto j = 0u; j < N; j++){
		for(auto i = j+1; i < N; i++){
			if(M[j][i] > thresh || M[j][i] < -thresh){
				join_proposal(i,j);
			}
		}
	}
}


/// Joins proposals parameters th1 and th2
void Update::join_proposal(unsigned int th1, unsigned int th2) 
{
	if(!model.param_vec[th1].prop_pos) return;
	if(!model.param_vec[th2].prop_pos) return;
	
	auto th1_p = UNSET, th2_p = UNSET; 
	for(auto p = 0u; p < proposal.size(); p++){
		const auto &prop = proposal[p];
		if(prop.type == PARAM_PROP){
			if(find_in(prop.param_list,th1) != UNSET) th1_p = p;
			if(find_in(prop.param_list,th2) != UNSET) th2_p = p;
		}
	}
	if(th1_p == UNSET || th2_p == UNSET) emsg("th unset");
				
	if(th1_p == th2_p) return;
	
	auto &prop1 = proposal[th1_p]; 
	auto &prop2 = proposal[th2_p];
	
	const auto &plist1 = prop1.param_list;
	const auto &plist2 = prop2.param_list;
	
	// Switches off existing proposals (if they have more than one variable)
	for(auto &prop : proposal){
		if(prop.type == PARAM_PROP || prop.type == MBP_PROP){
			auto si_lim = 1u; if(prop.type == MBP_PROP) si_lim = 0;
			if(plist1.size() > si_lim && equal_vec(prop.param_list,plist1)) prop.on = false;
			if(plist2.size() > si_lim && equal_vec(prop.param_list,plist2)) prop.on = false;
		}
	}
	
	auto param_list = combine(plist1,plist2);
	
	auto j_st = proposal.size();
	
	add_parameter_prop(param_list);
	
	for(auto j = j_st; j < proposal.size(); j++){
		auto &prop = proposal[j];
		cout << "JOIN ADD: " << prop.name << endl;
		prop.update_sampler(cor_matrix);
		
		prop.calculate_affect_spline(); 
	}
}


