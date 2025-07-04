// This file deals wiht the system state (i.e. the parameters and dynamics for all the species)

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "state.hh"
#include "utils.hh"
#include "matrix.hh"
#include "ind_ev_sampler.hh"

/// Initialises the model 
State::State(const Model &model) : model(model)
{
	dif_thresh = DIF_THRESH;
}


/// Initialises state;
void State::init()
{
	T = model.ntimepoint-1;
	
	timer.resize(TIMER_MAX,0);
	update_timer.resize(AFFECT_MAX,0);
	restore_timer.resize(AFFECT_MAX,0);

	dtimepoint.resize(T);
	for(auto ti = 0u; ti < T; ti++) dtimepoint[ti] = model.timepoint[ti+1] - model.timepoint[ti];
	
	spline_val.clear();
	spline_val.resize(model.spline.size());
	for(auto &spl : spline_val){
		spl.val.resize(T,UNSET);
		spl.store.resize(T,UNSET);
	}
	
	popnum_t.resize(T);

	nspecies = model.nspecies;
	for(auto p = 0u; p < nspecies; p++){
		vector <unsigned int> pop_affect;               // Lists populations affect by species
		for(auto po = 0u; po < model.pop.size(); po++){
			if(model.pop[po].sp_p == p) pop_affect.push_back(po);
		}
		
		StateSpecies ss(param_val,spline_val,model.eqn,model.param,model.param_vec,model.pop,model.species[p],model.genetic_data,model.details,model.timepoint,dtimepoint,pop_affect,model.mode,dif_thresh);
		species.push_back(ss);
	}
	
	initialise_update_ind_maps();
	initialise_local_map();
}


/// Simulate the state from parameter values, individual effects and initial conditions
void State::simulate(const vector <double> &param_value, const vector <InitCondValue> &initc_val) 
{
	param_val = param_value;

	spline_init(); 

	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		ssp.init_cond_val = initc_val[p];
		ssp.simulate_init();
		ssp.simulate_individual_init();
	}
		
	popnum_t[0] = model.calculate_popnum(species);  

	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		if(ssp.type == INDIVIDUAL){
			ssp.activate_initial_state(model.details.t_start,popnum_t);
		}
	}

	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		ssp.check(0,popnum_t);
	}

	simulate_iterate(0,T);
}


/// Simulate the state from parameter values, individual effects and initial conditions
void State::post_sim(const vector <double> &param_value, const Sample &samp)
{
	param_val = param_value;
	
	spline_init(); 
	
	auto t_start = model.details.ppc_t_start;
	
	auto t_end = model.details.ppc_t_end;
	auto ti_start = get_ti(t_start);
	auto ti_end = get_ti(t_end);
	
	for(auto p = 0u; p < model.species.size(); p++){
		auto &ssp = species[p];
		ssp.simulate_init();
		ssp.simulate_sample_init(ti_start,samp.species[p]);
	}
	
	popnum_t = model.calculate_popnum_t(species,ti_start+1);
	popnum_t.resize(T);
	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		switch(ssp.type){
		case INDIVIDUAL:
			ssp.activate_initial_state(t_start,popnum_t);
			break;
			
		case POPULATION:
			ssp.set_tnum_mean(ti_start,popnum_t);
			ssp.cpop = ssp.cpop_st[ti_start];
			ssp.cpop_st.pop_back();
			break;
		}
	}

	simulate_iterate(ti_start,ti_end);
}


/// Once a simulation is setup this iterates equations
void State::simulate_iterate(unsigned int ti_start, unsigned int ti_end) 
{
	for(auto ti = ti_start; ti < ti_end; ti++){	
		//if(true && ti%op_step == 0) print_cpop(ti);
	
		auto &pop = popnum_t[ti];
	
		pop = model.calculate_popnum(species);          // Calculates the population numbers
		
		auto pop_ind = calculate_pop_ind(); 
	
		for(auto p = 0u; p < nspecies; p++){
			auto &ssp = species[p];
			switch(ssp.type){
			case INDIVIDUAL:
				ssp.update_individual_based(ti,pop_ind,popnum_t);
				break;
				
			case POPULATION: 
				ssp.update_population_based(ti,model.details.stochastic,pop);
				break;
			}
			ssp.check(ti,popnum_t);
		}
	}

	if(model.mode == INF) ensure_all_ind_event();

	if(false){
		cout << "PRINT" << endl;	
		for(auto p = 0u; p < nspecies; p++){
			auto &ssp = species[p];
			for(auto i = 0u; i < ssp.individual.size(); i++){
			//for(auto i = 0u; i <10; i++){
				ssp.print_ind(i);
				ssp.print_event("start",ssp.individual[i]);
			}	
		}	
		emsg("done");
	}
	
	if(false){
		auto prob_sum = 0.0;
		for(const auto &ssp : species) prob_sum += ssp.prob_trans_tree;
		if(prob_sum != 0) cout << "Probability trans tree = " << prob_sum << endl;
		//emsg("do");
	}

	likelihood_from_scratch();
}


/// Ensures that all individuals have event
// (If the model has a source its possible simulation results in no event
void State::ensure_all_ind_event()
{
	for(auto p = 0u; p < species.size(); p++){
		const auto &sp = model.species[p];
		auto &ssp = species[p];
		
		if(sp.type == INDIVIDUAL && sp.contains_source){	
			IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,genetic_value.inf_node,param_val,spline_val,popnum_t);
	
			for(auto i = 0u; i < sp.nindividual_obs; i++){
				auto &ev = ssp.individual[i].ev;
				if(ev.size() == 0){
					const auto &indd = sp.individual[i];
					
					if(indd.obs.size() == 0) emsg("Should be non-zero");
					auto trange = ssp.source_time_range(indd);
					
					auto f = 0.9;
					auto t = (1-f)*trange.tmin+f*trange.tmax;
					
					auto c = (unsigned int)(ran()*sp.comp_gl.size());
					IndInfFrom inf_from;
					auto e_init = ssp.get_event(ENTER_EV,i,UNSET,UNSET,UNSET,c,t,inf_from);
			
					auto temp = 0.0;
					ev = ind_ev_samp.simulate_events(i,e_init,temp,indd.trig_ev_ref);
				}
			}
		}
	}
}


/// Calcualtes derived quantities
vector <DeriveOutput> State::derive_calculate() const
{
	vector <DeriveOutput> output;
	
	vector < vector < vector <double> > > der_value; 
	
	for(const auto &der : model.derive){
		DeriveOutput op;
		vector < vector <double> > d_value;
		if(der.func.on){
			auto value = calculate_df(der.func);
			auto val_str = compact_vector(value);
			op.value_str.push_back(val_str);
			d_value.push_back(value);
		}
		else{
			for(auto i = 0u; i < der.eq.size(); i++){
				auto &eqn = model.eqn[der.eq[i].eq_ref];
				vector <double> value;
				
				string val_str;
				if(der.time_dep == false){
					auto val = eqn.calculate_derive(UNSET,popnum_t,param_val,spline_val,der_value);
					val_str = tstr(val);
					value.push_back(val);
				}
				else{
					for(auto ti = 0u; ti < T; ti++){	
						auto val = eqn.calculate_derive(ti,popnum_t,param_val,spline_val,der_value);
						value.push_back(val);
					}
					val_str = compact_vector(value);
				}
				d_value.push_back(value);
				
				op.value_str.push_back(val_str);
			}
		}
		
		output.push_back(op);
		der_value.push_back(d_value);
	}
	
	return output;
}


/// Calculates derived functions
vector <double> State::calculate_df(const DerFunc &df) const
{
	vector <double> vec;
	vector < vector <double> > F, S, Sinv, V;
	
	auto &F_eq = df.F_eq;
	auto &S_eq = df.S_eq;
	auto &Feq_ref = df.Feq_ref;
	auto &Seq_ref = df.Seq_ref;
	
	auto N = F_eq.size();
	
	F.resize(N); S.resize(N); 
	for(auto c = 0u; c < N; c++){ F[c].resize(N); S[c].resize(N);}
	
	const auto &sp = model.species[0];
	const auto &ssp = species[0];
	 
	// Initial guess for eigen vector
	vector <double> evec(N);
	for(auto i = 0u; i < N; i++) evec[i] = ran();
		
	auto cpop = ssp.cpop_st;
	if(sp.type != POPULATION) cpop = ssp.ibm_cpop_st();
		
	auto type = df.type;
	
	// Computational approach to RN
	if(type == RNC || type == GTC){
		const auto &ref = df.ref;
		auto dt = model.details.dt;
		const auto &timepoint = model.timepoint;
		
		// Gets which individuals are in which populations
		auto pop_ind = calculate_pop_ind_total(); 

		// Gets periods over which individuals are infected
		auto inf_period = ssp.get_inf_period(ref);
		
		// Goes through infection events and allocates contributions to individuals
		for(const auto &ind : ssp.individual){
			for(const auto &ev : ind.ev){
				if(ev.type == M_TRANS_EV){
					const auto &tr = sp.tra_gl[ev.tr_gl];
					if(tr.i != UNSET && tr.f != UNSET){
						if(ref[tr.i] == UNSET && ref[tr.f] != UNSET){
							auto t = ev.t;
							auto ti = get_ti(t);
							
							const auto &eqn = model.eqn[tr.dist_param[0].eq_ref];
							const auto &lin = eqn.linearise;
							if(!lin.on) emsg("Linearisation should be on");
							
							// Generates a list of all individuals which could have cause infection
							vector <Poss> inf_pos;
							auto sum = 0.0;
							for(auto k = 0u; k < eqn.pop_ref.size(); k++){
								auto po = eqn.pop_ref[k];
								auto val = eqn.calculate_calculation(lin.pop_grad_calc[k],ti,param_val,spline_val);
								for(auto poss : pop_ind[ti][po]){
									poss.weight *= val;
									inf_pos.push_back(poss);
									sum += poss.weight;
								}
								
								for(auto &pos : inf_pos){
									auto i = pos.i;
									auto k = 0u;
									while(k < inf_period[i].size() && t > inf_period[i][k].end_inf) k++;	
									if(k == inf_period[i].size()) emsg("Could not find in period");
									
									auto val = pos.weight/sum;
									
									inf_period[i][k].num_inf += val;
									inf_period[i][k].num_inf_t += val*(t-inf_period[i][k].start);
								}					
							}
						}
					}
				}
			}						
		}
	
		vec.resize(T,0);
		vector <double> num(T,0);
		for(auto i = 0u; i < inf_period.size(); i++){
			for(const auto &ip : inf_period[i]){
				auto start = ip.start;
				auto end = ip.end;
				auto ww = 1.0;
				//auto ww = 1.0/(end-start); // This seems to give values too small
				auto val = ip.num_inf;
				if(type == GTC){
					val = ip.num_inf_t/(ip.num_inf+TINY);
				}				
				auto ti = get_ti(start);
				auto ti_end = get_ti(end);
				if(ti == ti_end){
					auto w = ww*(end-start);
					vec[ti] += val*w;
					num[ti] += w;
				}
				else{
					auto wi = ww*(timepoint[ti+1]-start);
					vec[ti] += val*wi;
					num[ti] += wi;
					
					for(auto tii = ti+1; tii < ti_end; tii++){
						auto w = ww*dt;
						vec[tii] += val*w;
						num[tii] += w;
					}
					
					auto wf = ww*(end-timepoint[ti_end]);
					vec[ti_end] += val*wf;
					num[ti_end] += wf;	
				}
			}
		}
		
		for(auto ti = 0u; ti < T; ti++) vec[ti] /= (num[ti]+TINY);
	
		if(false){
			for(auto i = 0u; i < inf_period.size(); i++){
				cout << i << ": ";
				for(auto inf : inf_period[i]){
					cout << inf.num_inf << "  ";
				}
				cout << endl;
			}
			emsg("do");
		}
	}
	else{
		for(auto ti = 0u; ti < T; ti++){	
			vector <double> Fvalue;
			for(auto i = 0u; i < Feq_ref.size(); i++){
				const auto &fe = Feq_ref[i];
				const auto &eqn = model.eqn[fe.e];
				auto val = eqn.calculate_calculation(eqn.linearise.pop_grad_calc[fe.i],ti,param_val,spline_val);
				Fvalue.push_back(val);
			}
			
			vector <double> Svalue;
			for(auto i = 0u; i < Seq_ref.size(); i++){
				const auto &se = Seq_ref[i];
				
				// Calculates any branching factor
				auto bp = 1.0;  
				if(se.m != UNSET){
					const auto &nmt = sp.nm_trans[se.m];
					const auto &eqn = model.eqn[nmt.bp_eq];
					
					bp = eqn.calculate_no_popnum(ti,param_val,spline_val);
					if(nmt.all_branches){
						auto bp_tot = 0.0;
						for(auto bp_eq : nmt.bp_all_eq){
							const auto &eqn = model.eqn[bp_eq];
							bp_tot += eqn.calculate_no_popnum(ti,param_val,spline_val);
						}
						bp /= bp_tot;
					}
				}
			
				const auto &eqn = model.eqn[se.e];
				auto val = eqn.calculate_no_popnum(ti,param_val,spline_val);
				
				switch(se.type){
				case EXP_RATE: case EXP_RATE_NM: 
					Svalue.push_back(bp*val);
					break;
				case EXP_MEAN: case EXP_MEAN_NM: case GAMMA: 
				case ERLANG: case LOG_NORMAL: case PERIOD: 
					Svalue.push_back(bp/(val+TINY));
					break;
					
				case WEIBULL: emsg("SHould not be weibull"); break; 
				}
			}
			
			for(auto c = 0u; c < N; c++){
				for(auto cc = 0u; cc < N; cc++){
					auto Fsum = 0.0;
					for(const auto &re : F_eq[c][cc]){
						if(type == RN || type == GT) Fsum += Fvalue[re.k]*cpop[0][re.c_from];
						else Fsum += Fvalue[re.k]*cpop[ti][re.c_from];
					}
					F[c][cc] = Fsum;
					
					auto Ssum = 0.0;
					for(const auto &re : S_eq[c][cc]){
						Ssum += Svalue[re.k]*re.sign;
					}
					S[c][cc] = Ssum;
				}
			}
			
			//print_matrix("F",F);
			//print_matrix("S",S);
			
			Sinv = invert_matrix(S);
			V = matrix_mult(F,Sinv);
			
			auto Rval = largest_eigenvalue(V,evec);
				
			switch(type){
			case RN: case RNE:
				vec.push_back(Rval);
				break;
				
			case GT: case GTE:
				{
					auto Z = matrix_mult(Sinv,evec);
					vector <double> tau(N,0);
					for(auto c = 0u; c < N; c++){
						for(auto e : df.tau_eq[c]){
							tau[c] += 1.0/Svalue[e];
						}
						tau[c] += 0.5/S[c][c];
					}
					
					auto sum_top = 0.0, sum_bot = 0.0;
					for(auto c = 0u; c < N; c++){
						for(auto cc = 0u; cc < N; cc++){
							sum_top += F[c][cc]*tau[cc]*Z[cc];
							sum_bot += F[c][cc]*Z[cc];
						}
					}
					vec.push_back(sum_top/(sum_bot+TINY));
				}
				break;
			
			default: emsg("op not here"); break;
			}
		}
	}
	vec.push_back(vec[vec.size()-1]);
	
	return vec;
}
	
	
/// Compactifies a vector into a string
string State::compact_vector(const vector <double> &value) const
{
	string st;
	auto nt = 0u;
	auto i = 0u;
	while(i < value.size()){
		auto val = value[i];
		auto ist = i;
		while(i+1 < value.size() && value[i+1] == val) i++;
		i++;
		
		if(ist != 0) st += "|";
		st += tstr(val);
		auto num = i-ist;
		if(num > 1) st += "*"+tstr(num);
		nt += num;
	}
	
	return st;
}


/// Works out values for the spline based on the parameter vector
void State::spline_init()
{
	for(auto j = 0u; j < model.spline.size(); j++){
		const auto &spline = model.spline[j];
		const auto &div = spline.div;
		auto &val = spline_val[j].val;
		
		if(spline.constant == false){
			for(auto ti = 0u; ti < T; ti++){
				auto f = div[ti].f;
				val[ti] = param_val[div[ti].th1]*f + param_val[div[ti].th2]*(1-f);
			}
		}
		else val = spline.const_val;
	}
}


/// Prints cpop to the console
void State::print_cpop(unsigned int ti) const 
{
	if(com_op == true) return;
	
	cout << model.timepoint[ti] << ": ";
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = model.species[p];
		const auto &ssp = species[p];
		if(nspecies > 1) cout << model.species[p].name << " - ";
		
		const auto C = sp.comp_gl.size();
		vector <double> num(C,0);
		auto n = 0u;
		for(auto c = 0u; c < C; c++){
			const auto &co = sp.comp_gl[c];
			if(co.erlang_hidden){ num[co.erlang_c_start] += ssp.cpop[c];}
			else{ num[c] += ssp.cpop[c]; n++;}
		}
		
		if(n <= 10){
			for(auto c = 0u; c < C; c++){
				const auto &co = sp.comp_gl[c];
				if(!co.erlang_hidden){
					cout << sp.comp_gl[c].name <<"=" << num[c] << " ";
				}
			}
		}
		else{
			// Looks at just the first classification
			auto cl = 0;
			
			const auto &claa = sp.cla[cl]; 
			const auto CC = claa.ncomp;
			
			auto nn = 0u;
			for(auto cc = 0u; cc < CC; cc++){
				if(!claa.comp[cc].erlang_hidden) nn++;
			}
			if(nn < 10){
				vector <double> numC(CC,0);
				for(auto c = 0u; c < C; c++){
					numC[sp.comp_gl[c].cla_comp[cl]] += num[c];
				}
				
				for(auto cc = 0u; cc < CC; cc++){
					if(!claa.comp[cc].erlang_hidden){
						cout << claa.comp[cc].name << "=" << numC[cc] << " ";
					}
				}
			}
			else{
				cout << "P=" << sum(num) << " ";
			}
		}
		
		cout << "   ";
	}
	cout << endl;
}


/// Calculates the likelihood for all parts of the model
void State::calculate_likelihood()
{
	prior_prob = model.prior_prob(param_val);
	
	dist_prob = model.dist_prob(param_val);
	
	spline_prior = model.spline_prior(param_val);
	
	auto all_time = seq_vec(T);
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &sp = model.species[p];
		auto &ssp = species[p];

		ssp.calculate_obs_eqn(seq_vec(sp.obs_eqn.size()));
		
		switch(ssp.type){
		case INDIVIDUAL:
			{
				for(auto m = 0u; m < sp.obs_trans_eqn.size(); m++){
					auto temp = 0.0;
					ssp.likelihood_unobs_trans(m,all_time,temp);
				}
					
				for(auto g = 0u; g < ssp.ind_eff_group_sampler.size(); g++){
					auto temp = 0.0;
					ssp.likelihood_indeff_group(g,temp);
				}
			
				ssp.Li_markov.resize(ssp.N);
				for(auto e = 0u; e < ssp.N; e++){
					auto &me_vari = ssp.markov_eqn_vari[e];
					ssp.Li_markov[e].resize(me_vari.div.size());
					auto temp = 0.0;
					auto list = seq_vec(me_vari.div.size());
					ssp.likelihood_markov(e,list,temp);
				}
			
			
				for(auto i = 0u; i < sp.nm_trans.size(); i++){
					auto temp = 0.0;
					ssp.likelihood_nm_trans(i,all_time,popnum_t,temp);
					ssp.likelihood_nm_trans_bp(i,all_time,popnum_t,temp);
				}
				
				for(auto i = 0u; i < sp.nm_trans_incomp.size(); i++){
					auto temp = 0.0;
					ssp.likelihood_nm_trans_incomp(i,all_time,popnum_t,temp);
				}
				
				ssp.Li_obs_ind.resize(sp.nindividual_in);	
				auto temp = 0.0;
				ssp.likelihood_obs_ind(seq_vec(sp.nindividual_in),temp);
			}
			break;
			
		case POPULATION:
			{	
				if(model.details.stochastic == true) ssp.likelihood_pop(popnum_t);
			}
			break;
		}
		
		ssp.calculate_pop_data_cgl_trgl();
		
		auto temp = 0.0;
		ssp.Li_obs_pop.resize(sp.pop_data.size());
		ssp.likelihood_obs_pop(seq_vec(sp.pop_data.size()),temp);
		
		ssp.Li_obs_pop_trans.resize(sp.pop_trans_data.size());
		ssp.likelihood_obs_pop_trans(seq_vec(sp.pop_trans_data.size()),temp);
		
		ssp.likelihood_init_cond(temp);
	}
	
	calculate_like();
}


/// Calculates the overall likelihood 
void State::calculate_like()
{
	like.prior = sum(prior_prob);
	like.dist = sum(dist_prob);
	like.spline_prior = sum(spline_prior);
	
	like.init_cond = 0;
	like.init_cond_prior = 0;
	prior_init_cond(like.init_cond_prior);
	like.markov = 0;
	like.nm_trans = 0;
	like.genetic_process = likelihood_genetic_process();
	like.genetic_obs = likelihood_genetic_obs(genetic_value.gen_dif);
	like.ie = 0;
	like.obs = 0;
	for(auto p = 0u; p < nspecies; p++){
		auto &ssp = species[p];
		
		like.obs += sum(ssp.Li_obs_ind) + sum(ssp.Li_obs_pop) + sum(ssp.Li_obs_pop_trans)
						+ sum_mult(ssp.obs_trans_eqn_value,ssp.obs_trans_eqn_num);
	
		like.init_cond += ssp.Li_init_cond;
	
		switch(ssp.type){
		case INDIVIDUAL:
			like.markov += sum(ssp.Li_markov);
			
			{
				auto Li_nm = 0.0; 
				
				for(const auto &ind : ssp.individual){
					for(const auto &ev : ind.ev){
						if(ev.type == NM_TRANS_EV){
							Li_nm += ev.Li + ev.Li_bp;
						}
					}
					
					for(const auto &inmtr : ind.incomp_ref){
						if(inmtr.on) Li_nm += inmtr.Li;
					}
				}
				
				like.nm_trans += Li_nm;
			}
			
			like.ie += sum(ssp.Li_ie);
			break;
			
		case POPULATION:
			like.markov += sum(ssp.Li_markov_pop);
			break;
		}
	}
}


/// Adds a change in likelihood to the state
void State::add_like(Like like_ch)
{
	like.init_cond += like_ch.init_cond;
	like.init_cond_prior += like_ch.init_cond_prior;
	like.obs += like_ch.obs;
	like.prior += like_ch.prior;
	like.dist += like_ch.dist;
	like.ie += like_ch.ie;
	like.spline_prior += like_ch.spline_prior;
	like.markov += like_ch.markov;
	like.nm_trans += like_ch.nm_trans;
	like.genetic_process += like_ch.genetic_process;
	like.genetic_obs += like_ch.genetic_obs;
}


/// Updates splines
void State::update_spline(const vector <AffectLike> &affect_spline)
{
	vector <double> store;
	for(const auto &aspline : affect_spline){
		if(aspline.type != SPLINE_AFFECT) emsg("SHould be spline affect");
	
		auto j = aspline.num; 
		
		const auto &spline = model.spline[j];
		const auto &div = spline.div;
			
		if(spline.constant == true) emsg("SHould not recalculate");
		
		auto &val = spline_val[j].val;
		auto &store = spline_val[j].store;
		
		for(auto ti : aspline.list){
			auto f = div[ti].f;
			store[ti] = val[ti];
			val[ti] = param_val[div[ti].th1]*f + param_val[div[ti].th2]*(1-f);
		}
	}
}


/// Restores splines to their original value
void State::restore_spline(const vector <AffectLike> &affect_spline)
{
	for(const auto &aspline : affect_spline){
		if(aspline.type != SPLINE_AFFECT) emsg("SHould be spline affect");
	
		auto j = aspline.num; 
		auto &val = spline_val[j].val;
		auto &store = spline_val[j].store;
		for(auto ti : aspline.list){ val[ti] = store[ti]; store[ti] = UNSET;}
	}
}


/// Restores splines to their original value
void State::remove_store_spline(const vector <AffectLike> &affect_spline)
{
	for(const auto &aspline : affect_spline){
		if(aspline.type != SPLINE_AFFECT) emsg("SHould be spline affect");
	
		auto j = aspline.num; 
		auto &store = spline_val[j].store;
		for(auto ti : aspline.list){ store[ti] = UNSET;}
	}
}


/// Updates state based on a change in parameters
Like State::update_param(const vector <AffectLike> &affect_like, const vector <double> &param_store)
{	
	Like like_ch;
	
	change.clear();
	
	auto K = affect_like.size();
	for(auto i = 0u; i < K; i++){
		const auto &alike = affect_like[i];
		
		update_timer[alike.type] -= clock();
		
		switch(alike.type){
		case OBS_EQN_AFFECT:
			{
				change_add(species[alike.num].calculate_obs_eqn(alike.list));
			}
			break;
			
		case LIKE_UNOBS_TRANS_AFFECT:
			{
				change_add(species[alike.num].likelihood_unobs_trans(alike.num2,alike.list,like_ch.obs));
			}
			break;

		case POP_DATA_CGL_TGL_AFFECT:
			{
				change_add(species[alike.num].calculate_pop_data_cgl_trgl());
			}
			break;
		
		case LIKE_INIT_COND_AFFECT:
			{
				change_add(species[alike.num].likelihood_init_cond(like_ch.init_cond));
			}
			break;
			
		case PRIOR_INIT_COND_AFFECT:
			{
				change_add(prior_init_cond(like_ch.init_cond_prior));
			}
			break;
			
		case LIKE_OBS_IND_AFFECT:
			{
				change_add(species[alike.num].likelihood_obs_ind(alike.list,like_ch.obs));
			}
			break;
			
		case LIKE_OBS_POP_AFFECT:
			{
				change_add(species[alike.num].likelihood_obs_pop(alike.list,like_ch.obs));
			}
			break;
			
		case LIKE_OBS_POP_TRANS_AFFECT:
			{
				change_add(species[alike.num].likelihood_obs_pop_trans(alike.list,like_ch.obs));
			}
			break;
		
		case MARKOV_POP_AFFECT:
			{
				change_add(species[alike.num].likelihood_pop_change(alike.num2,alike.list,popnum_t,like_ch.markov));
			}
			break;
			
		case LIKE_IE_AFFECT:
			{
				change_add(species[alike.num].likelihood_indeff_group(alike.num2,like_ch.ie));
			}
			break;
			
		case EXP_IE_AFFECT:
			{
				change_add(species[alike.num].recalculate_exp_ie(alike.num2));
			}
			break;
			
		case OMEGA_AFFECT:
			{
				change_add(species[alike.num].calculate_omega(alike.num2));
			}
			break;
			
		case NM_TRANS_AFFECT:
			{
				change_add(species[alike.num].likelihood_nm_trans(alike.num2,alike.list,popnum_t,like_ch.nm_trans));
			}
			break;
			
		case NM_TRANS_BP_AFFECT:
			{
				change_add(species[alike.num].likelihood_nm_trans_bp(alike.num2,alike.list,popnum_t,like_ch.nm_trans));
			}
			break;
			
		case NM_TRANS_INCOMP_AFFECT:
			{
				change_add(species[alike.num].likelihood_nm_trans_incomp(alike.num2,alike.list,popnum_t,like_ch.nm_trans));
			}
			break;
			
		case POP_AFFECT:
			{	
				change_add(model.recalculate_population(popnum_t,alike.list,species));
			}
			break;
			
		case SPLINE_PRIOR_AFFECT:
			{
				change_add(model.recalculate_spline_prior(alike.num,spline_prior,param_val,like_ch.spline_prior));
			}
			break;
			
		case PRIOR_AFFECT:
			{
				change_add(model.recalculate_prior(alike.num,prior_prob,param_val,like_ch.prior));
			}
			break;
			
		case DIST_AFFECT:
			{
				change_add(model.recalculate_prior(alike.num,dist_prob,param_val,like_ch.dist));
			}
			break;
		
		case SPLINE_AFFECT: emsg("SHould not be spline affect"); break;
			
		case DIV_VALUE_AFFECT:   // Updates div.value on Markov transition
			{
				auto p = alike.num, e = alike.num2;
				change_add(species[p].likelihood_markov_value(e,alike.list,popnum_t));
			} 
			break;
			
		case DIV_VALUE_FAST_AFFECT:   // Updates div.value on Markov transition
			{
				auto p = alike.num;
				change_add(species[p].likelihood_markov_value_fast(alike.me_list,alike.list,popnum_t,param_store,spline_val));
			}
			break;
			
		case DIV_VALUE_LINEAR_AFFECT:   // Updates div.value on Markov transition
			{
				auto p = alike.num;
				change_add(species[p].likelihood_markov_value_linear(alike.list,alike.linear_prop,popnum_t));
			}
			break;
		
		case INDFAC_INT_AFFECT:      // Updates div.indfac_int on Markov transition
			{
				auto p = alike.num;
				change_add(species[p].likelihood_indfac_int());
			}
			break;
			
		case EXP_FE_AFFECT:      // Updates exp_fe on individual
			{
				auto p = alike.num, f = alike.num2;
				change_add(species[p].set_exp_fe(f));
			}
			break;
			
		case MARKOV_LIKE_AFFECT: // Updates the Markov liklihood
			{
				auto p = alike.num, e = alike.num2;	
				change_add(species[p].likelihood_markov(e,alike.list,like_ch.markov));
			}
			break;

		case LIKE_GENETIC_PROCESS_AFFECT:
			{
				like_ch.genetic_process = likelihood_genetic_process() - like.genetic_process;
				change_add();
			}
			break;
			
		case GENETIC_VALUE_AFFECT:
			{
				change_add(set_genetic_param());
			}
			break;
			
		case LIKE_GENETIC_OBS_AFFECT:
			{
				like_ch.genetic_obs = likelihood_genetic_obs(genetic_value.gen_dif) - like.genetic_obs;
				change_add();
			}
			break;
			
		case IIF_W_AFFECT:
			{	
				change_add(recalculate_iif_w());
			}
			break;
			
		case POPNUM_IND_W_AFFECT:
			{	
				change_add(recalculate_popnum_ind_w(alike.num));
			}
			break;
			
		case AFFECT_MAX: break;
		}
		
		update_timer[alike.type] += clock();
	}
	
	return like_ch;
}


/// Restores state to as before proposal
void State::restore(const vector <AffectLike> &affect_like)
{	
	for(auto i = 0u; i < change.size(); i++){
		const auto &alike = affect_like[i];
		const auto &vec = change[i].vec;
		const auto &li = alike.list;
		
		restore_timer[alike.type] -= clock();
		
		switch(alike.type){	
		case OBS_EQN_AFFECT:
			{
				auto p = alike.num;
				for(auto k = 0u; k < li.size(); k++) species[p].obs_eqn_value[li[k]] = vec[k];	
			}
			break;
			
		case LIKE_UNOBS_TRANS_AFFECT:
			{
				auto p = alike.num;
				auto &otev = species[p].obs_trans_eqn_value[alike.num2];
				for(auto k = 0u; k < li.size(); k++) otev[li[k]] = vec[k];	
			}
			break;
			
		case POP_DATA_CGL_TGL_AFFECT:
			{
				species[alike.num].restore_pop_data_cgl_trgl(vec);
			}
			break;
		
		case LIKE_INIT_COND_AFFECT:
			{
				species[alike.num].Li_init_cond = vec[0];
			}
			break;
			
		case PRIOR_INIT_COND_AFFECT:
			{
				like.init_cond_prior = vec[0]; 
			}
			break;
			
		case LIKE_OBS_IND_AFFECT:
			{
				auto p = alike.num;
				for(auto k = 0u; k < li.size(); k++) species[p].Li_obs_ind[li[k]] = vec[k];	
			}
			break;
			
		case LIKE_OBS_POP_AFFECT:
			{
				auto p = alike.num;
				auto k = 0u;
				for(auto i : li){ 
					species[p].Li_obs_pop[i] = vec[k]; k++;
					species[p].pop_data_num[i] = vec[k]; k++;
				}
				if(k != vec.size()) emsg("Problem with k");
			}
			break;
			
		case LIKE_OBS_POP_TRANS_AFFECT:
			{
				auto p = alike.num;
				auto k = 0u;
				for(auto i : li){ 
					species[p].Li_obs_pop_trans[i] = vec[k]; k++;
					species[p].pop_trans_data_num[i] = vec[k]; k++;
				}
				if(k != vec.size()) emsg("Problem with k");
			}
			break;
		
		case MARKOV_POP_AFFECT:
			{
				species[alike.num].likelihood_pop_change_restore(alike.num2,alike.list,vec);
			}
			break;
			
		case LIKE_IE_AFFECT:
			{
				species[alike.num].Li_ie[alike.num2] = change[i].num;
			}
			break;
			
		case EXP_IE_AFFECT:
			{
				species[alike.num].recalculate_exp_ie_restore(alike.num2,vec);
			}
			break;
			
		case OMEGA_AFFECT:
			{
				species[alike.num].calculate_omega_restore(alike.num2,vec);
			}
			break;
			
		case NM_TRANS_AFFECT:
			{	
				species[alike.num].likelihood_nm_trans_restore(alike.num2,li,vec);				
			}
			break;
			
		case NM_TRANS_BP_AFFECT:
			{
				species[alike.num].likelihood_nm_trans_bp_restore(alike.num2,li,vec);					
			}
			break;
			
		case NM_TRANS_INCOMP_AFFECT:
			{	
				species[alike.num].likelihood_nm_trans_incomp_restore(alike.num2,li,vec);				
			}
			break;
			
		case POP_AFFECT:
			{
				model.recalculate_population_restore(popnum_t,alike.list,vec);
			}
			break;
			
		case SPLINE_PRIOR_AFFECT:
			{
				spline_prior[alike.num] = change[i].num;
			}
			break;
			
		case PRIOR_AFFECT:
			{
				prior_prob[alike.num] = change[i].num;
			}
			break;
			
		case DIST_AFFECT:
			{
				dist_prob[alike.num] = change[i].num;
			}
			break;
			
		case SPLINE_AFFECT: emsg("Should not be spline affect"); break;
		
		case DIV_VALUE_AFFECT:                 // Restore div.value on Markov transition
			{
				auto p = alike.num, e = alike.num2;
				auto &div = species[p].markov_eqn_vari[e].div;
				for(auto k = 0u; k < li.size(); k++) div[li[k]].value = vec[k];
			}
			break;
			
		case DIV_VALUE_FAST_AFFECT:  
			{
				auto &mev = species[alike.num].markov_eqn_vari;
				auto m = 0u;
				for(auto k = 0u; k < li.size(); k++){
					for(auto e : alike.me_list){
						mev[e].div[li[k]].value = vec[m];
						m++;
					}
				}
			}
			break;
			
		case DIV_VALUE_LINEAR_AFFECT:   // Restore div.value on Markov transition
			{
				auto p = alike.num;
				species[p].likelihood_markov_value_linear_restore(alike.list,alike.linear_prop,vec);
			}
			break;
			
		case INDFAC_INT_AFFECT:      // Restore div.indfac_int on Markov transition
			{
				auto p = alike.num;
				species[p].likelihood_indfac_int_restore(vec);
			}
			break;
			
		case EXP_FE_AFFECT:  
			{
				auto p = alike.num, f = alike.num2;
				
				species[p].set_exp_fe_restore(f,vec);
			}
			break;
			
		case MARKOV_LIKE_AFFECT:     // Updates the Markov liklihood
			{
				auto p = alike.num, e = alike.num2;
				for(auto k = 0u; k < li.size(); k++) species[p].Li_markov[e][li[k]] = vec[k];	
			}
			break;
			
		case LIKE_GENETIC_PROCESS_AFFECT:	case LIKE_GENETIC_OBS_AFFECT:
			break;
			
		case GENETIC_VALUE_AFFECT:
			{
				genetic_value.mut_rate = vec[0];
				genetic_value.seq_var = vec[1];
			}
			break;
			
			
		case IIF_W_AFFECT:
			recalculate_iif_w_undo(vec);
			break;
			
		case POPNUM_IND_W_AFFECT:
			recalculate_popnum_ind_w_undo(alike.num,vec);
			break;
			
		case AFFECT_MAX: break;
		}
		
		restore_timer[alike.type] += clock();
	}
}


/// Adds a change to the state
void State::change_add(const vector <double> &vec)
{
	Change ch; ch.vec = vec; ch.num = UNSET;
	change.push_back(ch);
}

			
/// Adds a change to the state
void State::change_add(double num)
{
	Change ch; ch.num = num;
	change.push_back(ch);
}


/// Adds a change to the state
void State::change_add()
{
	Change ch; ch.num = UNSET;
	change.push_back(ch);
}


/// Sets up everything needed to represent the Markov likelihood
void State::likelihood_from_scratch()
{
	// Finds the populations
	popnum_t = model.calculate_popnum_t(species);

	// Sets up non-markovian transitions
	for(auto p = 0u; p < nspecies; p++){ 
		auto &ssp = species[p];
		ssp.setup_nm_trans(popnum_t);
		ssp.calculate_obs_trans_eqn_num();
	}
	
	calculate_popnum_ind();
	
	// Initialises markov_eqn with value calculated and inffac_int set to zero
	for(auto p = 0u; p < species.size(); p++){
		auto &sp = model.species[p];
		auto &ssp = species[p];
		
		switch(sp.type){
		case INDIVIDUAL:
			ssp.calculate_N_unobs();
			
			ssp.set_event_observed();
	
			for(auto g = 0u; g < sp.ind_eff_group.size(); g++){
				ssp.calculate_omega(g);
			}
			
			ssp.likelihood_indfac_int();
			ssp.likelihood_ind_trans();
			for(auto e = 0u; e < ssp.N; e++){
				auto list = seq_vec(ssp.markov_eqn_vari[e].div.size());
				ssp.likelihood_markov_value(e,list,popnum_t);
			}
			break;
			
		case POPULATION:
			ssp.Li_markov_pop.clear();
			ssp.Li_markov_pop.resize(sp.tra_gl.size());
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++) ssp.Li_markov_pop[tr].resize(T,0);
			break;
		}
		
		ssp.Li_init_cond = 0.0;
	}

	if(model.trans_tree) setup_transtree();
	
	if(model.genetic_data.on){
		set_genetic_param();
		sample_genetic_value();
	}
	
	calculate_likelihood();
}


/// Resamples individual using the observation sampler (this gets fixed events correct)
void State::resample_ind(bool do_pl)
{
	print_diag("Start resample...");
	
	auto popnum_t = model.calculate_popnum_t(species);
	
	auto pl = false; 
	if(do_pl) pl = true;

	if(pl) check(" before resample");
			
	for(auto p = 0u; p < species.size(); p++){
		auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			auto &ssp = species[p];
	
			IndEvSampler ind_ev_samp(ssp.markov_eqn_vari,ssp.individual,model.details,sp,ssp.obs_eqn_value,ssp.obs_trans_eqn_value,model.eqn,genetic_value.inf_node,param_val,spline_val,popnum_t);
			
			ind_ev_samp.setup_nm();
		
			for(auto i = 0u; i < sp.nindividual_in; i++){
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(ind_ev_samp.needed(i,cl) == true){
						if(pl){
							cout << endl << endl << endl << endl;
							cout << ssp.individual[i].name << endl;
							cout << "Before:" << endl; ssp.print_event(ssp.individual[i].ev);
						}
						
						ind_ev_samp.generate_ind_obs_timeline();
						auto loop_max = 1u;
						for(auto loop = 0u; loop < loop_max; loop++){
							auto probif = 0.0;
							auto ev_new = ind_ev_samp.sample_events(probif);
				
							if(pl){
								cout << endl << endl << endl << endl;
								cout << ssp.individual[i].name << endl;
								cout << "Before:" << endl; ssp.print_event(ssp.individual[i].ev);
								cout << endl;
								cout << "Propose:" << endl; ssp.print_event(ev_new);
								
								cout << "Illegal:" << ind_ev_samp.illegal << endl;
							}
						
							if(ind_ev_samp.illegal == false){
								auto gc = update_tree(p,i,ev_new);
								if(gc.type != GENCHA_FAIL){
									auto like_ch = update_ind(p,i,ev_new,UP_SINGLE);
									if(pl){
										cout << endl;
										cout << "After:" << endl; ssp.print_event(ssp.individual[i].ev);
									}
							
									add_like(like_ch);
									gen_change_update(gc); 	
									if(sp.trans_tree) update_popnum_ind(p,i);
									break;
								}
							}
						}
					}
				}
			}
		}
	}	
}

 
/// Generates a particle from the state
Particle State::generate_particle(unsigned int s, unsigned int chain, bool store_state) const
{
	Particle part;
	part.s = s; part.chain = chain;
	part.param_val_prop = model.get_param_val_prop(param_val);
	part.like = like;
	part.dir_out = derive_calculate();
	
	for(auto p = 0u; p < species.size(); p++){
		const auto &ssp = species[p];
		ParticleSpecies part_sp;
		part_sp.init_cond_val = ssp.init_cond_val;
		
		if(store_state){
			part_sp.trans_num = ssp.trans_num;
			part_sp.individual = ssp.individual;
			if(cum_diag && model.mode == INF) ssp.calc_trans_diag(part_sp,popnum_t);
		}
		part_sp.nindividual = ssp.individual.size();
		
		part.species.push_back(part_sp);
	}
		
	if(model.trans_tree){ 
		part.inf_origin = genetic_value.inf_origin;
		part.inf_node = genetic_value.inf_node;
	}
	
	if(model.trans_tree){                            // Outputs trans_tree statistics
		auto &tts = part.trans_tree_stats;
		const auto &gv = genetic_value;
		tts.N_mut_tree = 0u;
		tts.N_mut_origin = 0u;
		
		if(model.genetic_data.on){
			for(const auto &in : gv.inf_node){
				for(const auto &iev : in.inf_ev){
					tts.N_mut_tree += iev.mut_num;
				}
			}
			
			for(const auto &o : gv.inf_origin){
				tts.N_mut_origin += o.mut_num;
			}
		}
		
		tts.t_root = LARGE;
		for(const auto &in : gv.inf_node){
			if(in.t_start < tts.t_root) tts.t_root = in.t_start;
		}	
		
		tts.N_origin = gv.inf_origin.size();
		tts.N_inf = gv.inf_node.size();
		tts.N_unobs = gv.nobs_not_infected;
	}

	return part;
}


/// Generates a state from a particle 
void State::set_particle(const Particle &part, bool calc_like)
{
	param_val = model.get_param_val(part.param_val_prop);
	like = part.like;
	
	spline_init();

	for(auto p = 0u; p < species.size(); p++){
		auto &ssp = species[p];
		const auto &part_sp = part.species[p];
		ssp.init_cond_val = part_sp.init_cond_val;
		ssp.trans_num = part_sp.trans_num;
		ssp.individual = part_sp.individual;
		if(ssp.type == POPULATION) ssp.set_cpop_st();
	}
	
	if(model.genetic_data.on){
		genetic_value.inf_origin = part.inf_origin;
		genetic_value.inf_node = part.inf_node;
	}
	
	if(calc_like) likelihood_from_scratch();
	else popnum_t = model.calculate_popnum_t(species);
}


/// Calculate the prior probability of the initial condition
vector <double> State::prior_init_cond(double &like_ch)
{
	vector <double> store;

	if(model.mode == SIM) return store;
	
	store.push_back(like.init_cond_prior);
	
	auto prior = 0.0;
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		auto &ssp = species[p];
		
		if(sp.init_cond.type == INIT_POP_DIST){	
			const auto &ic = sp.init_cond;
			const auto &icv = ssp.init_cond_val;
	
			auto foc_cl = ic.focal_cl;

			if(foc_cl == UNSET){
				prior += prior_probability(icv.N_total,ic.pop_prior,param_val,model.eqn);
	
				prior += dirichlet_probability(icv.frac,ic.alpha);
			}
			else{                                        // A focal classification is set
				const auto &claa = sp.cla[foc_cl];
				for(auto c = 0u; c < claa.ncomp; c++){
					prior += prior_probability(icv.N_focal[c],ic.comp_prior[c],param_val,model.eqn);
				}
				
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(cl != foc_cl){
						prior += dirichlet_probability(icv.frac_focal[cl],ic.alpha_focal[cl]);
					}
				}
			}
		}
	}
	
	like_ch += prior-like.init_cond_prior;
	
	return store;
}


/// Updates the individual sampler (used to sample either the compartment of entry or insertion of source events
void State::update_individual_sampler()
{
	timer[UPDATE_SAMPLER_TIMER] -= clock();
	for(auto p = 0u; p < model.nspecies; p++){ 
		const auto &sp = model.species[p];
		if(sp.type == INDIVIDUAL){
			auto &ssp = species[p];
			ssp.source_sampler.update(sp.nindividual_obs,ssp.individual,ssp.markov_eqn_vari,sp.contains_source);                            // Optimises individual sampler
			
			ssp.update_rate_mean(popnum_t);
		}
	}
	
	timer[UPDATE_SAMPLER_TIMER] += clock();
}


/// Updates populations numbers in the initial vlaue popnum_t based on a change in the initial conditions
void State::update_popnum_t_init(unsigned int p, vector <unsigned int> cnum_i, vector <unsigned int> cnum_f)
{
	const auto &sp = model.species[p];
	
	for(auto c = 0u; c < sp.comp_gl.size(); c++){
		double d = double(cnum_f[c])-double(cnum_i[c]);
		if(d != 0){
			for(const auto &pr : sp.comp_gl[c].pop_ref){
				auto po = pr.po;
				const auto &te = model.pop[po].term[pr.index]; if(te.c != c) emsg("c wrong");
				popnum_t[0][po] += d*te.w;
			}
		}
	}
}


/// Initialises a way to restore system back to before individual changes
void State::back_init()
{
	back_pop.clear();
	for(auto &ssp : species){
		ssp.back.clear();
		ssp.back_ind.clear();
	}
}


/// Restores the state to the orignal
void State::restore_back()
{
	for(auto &ssp : species) ssp.restore_back();
	
	for(const auto &bp : back_pop){
		switch(bp.type){
		case POP_NUM_T:
			{
				auto ti = bp.i, ti_next = bp.j, k = bp.k;
				auto val = bp.value;
				for(auto tii = ti; tii < ti_next; tii++) popnum_t[tii][k] -= val;
			}
			break;
		}
	}
}


/// Gets ti div time from an actual time
unsigned int State::get_ti(double t) const
{
	return (unsigned int)(OVER_ONE*(t-model.details.t_start)/model.details.dt);
}


/// Gets param_val_prop from a full parameter vector (i.e. removes reparam)
vector <double> State::get_param_val_prop() const
{
	return model.get_param_val_prop(param_val);
}


/// Gets the estimated rate for a 
double State::get_trans_rate_est(unsigned int p, unsigned int tr, unsigned int ti) const
{
	const auto &tra = model.species[p].tra_gl[tr];
	const auto &eq = model.eqn[tra.dist_param[0].eq_ref];
		
	switch(tra.type){
	case EXP_RATE: case EXP_RATE_NM:
		return eq.calculate(ti,popnum_t[ti],param_val,spline_val); 
			
	case EXP_MEAN: case EXP_MEAN_NM: case GAMMA: case ERLANG: case LOG_NORMAL: 
	case PERIOD: case WEIBULL:
		return 1.0/(TINY+eq.calculate(ti,popnum_t[ti],param_val,spline_val)); 
	}
	return UNSET;
}
	
			
/// Gets transition rates  
vector < vector <double> > State::get_population_rates(unsigned int p) const
{
	const auto &sp = model.species[p];
	
	vector < vector <double> > rate;
	rate.resize(sp.tra_gl.size());
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		rate[tr].resize(T);
		
		const auto &tra = sp.tra_gl[tr];
		const auto &eq = model.eqn[tra.dist_param[0].eq_ref];
		if(sp.type == INDIVIDUAL && eq.pop_ref.size() == 0){
			for(auto ti = 0u; ti < T; ti++){		
				rate[tr][ti] = 0;
			}
		}
		else{
			if(eq.time_vari){
				for(auto ti = 0u; ti < T; ti++){		
					rate[tr][ti] = get_trans_rate_est(p,tr,ti);
				}
			}
			else{
				auto va = get_trans_rate_est(p,tr,0);
				for(auto ti = 0u; ti < T; ti++){		
					rate[tr][ti] = va;
				}
			}
		}
	}
	
	return rate;
}


/// Calculates the individuals which are associated with a population
vector < vector <Poss> > State::calculate_pop_ind() const 
{
	vector < vector <Poss> > pop_ind;
	
	auto p = 0u; while(p < nspecies && model.species[p].trans_tree == false) p++;
	if(p == nspecies) return pop_ind;
	
	pop_ind.resize(model.pop.size());
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &ssp = species[p];
		if(ssp.type == INDIVIDUAL){
			const auto &sp = model.species[p];
			for(auto i = 0u; i < ssp.individual.size(); i++){
				const auto &ind = ssp.individual[i];
				
				auto c = ssp.ind_sim_c[i];
				if(c != UNSET){
					for(const auto &pr : sp.comp_gl[c].pop_ref){
						const auto &po = model.pop[pr.po];	
						if(po.term[pr.index].c != c) emsg("Problem");
					
						auto num = 1.0; 
						if(po.ind_variation){
							for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
							for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
						}
				
						Poss poss; poss.i = i; poss.weight = po.term[pr.index].w*num;
						
						pop_ind[pr.po].push_back(poss);
					}
				}
			}
		}
	}
	
	return pop_ind;
}


/// Calculates the individuals which are associated with a population for all time
vector < vector < vector <Poss> > > State::calculate_pop_ind_total() const 
{
	vector < vector < vector <Poss> > > pop_ind;
	
	const auto &timepoint = model.timepoint;
	
	pop_ind.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		pop_ind[ti].resize(model.pop.size());
	}
	
	for(auto p = 0u; p < nspecies; p++){
		const auto &ssp = species[p];
		if(ssp.type == INDIVIDUAL){
			const auto &sp = model.species[p];
			for(auto i = 0u; i < ssp.individual.size(); i++){
				const auto &ind = ssp.individual[i];
				
				auto c = UNSET;
				auto ti = 0u;
				for(auto e = 0u; e <= ind.ev.size(); e++){
					double t;
					if(e < ind.ev.size()) t = ind.ev[e].t;
					else t = timepoint[T];
					
					if(c != UNSET){
						auto ti_start = ti;
						for(const auto &pr : sp.comp_gl[c].pop_ref){
							const auto &po = model.pop[pr.po];	
							if(po.term[pr.index].c != c) emsg("Problem");
						
							auto num = 1.0; 
							if(po.ind_variation){
								for(auto ie : po.ind_eff_mult) num *= ind.exp_ie[ie];
								for(auto fe : po.fix_eff_mult) num *= ind.exp_fe[fe];
							}
					
							Poss poss; poss.i = i; poss.weight = po.term[pr.index].w*num;
							
							ti = ti_start;
							while(timepoint[ti] < t){
								pop_ind[ti][pr.po].push_back(poss);
								ti++;
							}
						}
					}
					while(timepoint[ti] < t){
						ti++;
					}
					
					if(e < ind.ev.size()) c = ind.ev[e].c_after;
				}
			}
		}
	}
	
	if(false){
		for(auto ti = 0u; ti < T; ti++){
			for(const auto &pi : pop_ind[ti]){
				cout << timepoint[ti] << " " << pi.size() << " pop" << endl;
				
				for(auto poss : pi){
					cout << poss.i << " "<< poss.weight << "  iw" << endl;
				}
			}
		}
		emsg("do");
	}
	
	return pop_ind;
}

