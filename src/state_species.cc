// This deals with general state species functions and population-based updates

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
// NOTES
// markov_eqn list all transitions which are Markovian. Attached to this list are the individuals 

using namespace std;

#include "state_species.hh"
#include "utils.hh"
#include "matrix.hh"


/// Initialises the model model.param,model.param_vec
StateSpecies::StateSpecies(const vector <double> &param_val, const vector <SplineValue> &spline_val, const vector <Equation> &eqn, const vector <Param> &param, const vector <ParamVecEle> &param_vec, const Species &sp, const Details &details, const	vector <double> &timepoint, const	vector <double> &dtimepoint, const vector <unsigned int> &pop_affect_, Operation mode_) : param_val(param_val), spline_val(spline_val), eqn(eqn), param(param), param_vec(param_vec), sp(sp), details(details), timepoint(timepoint), dtimepoint(dtimepoint) 
{
	pop_affect = pop_affect_;
	
	N = sp.markov_eqn.size();
	T = timepoint.size()-1;
	ti_sim = 0;
	
	mode = mode_;
	
	nnode = sp.markov_tree.node.size();
	
	initialise_arrays();
}


/// Initialises all the arrays which will be used
void StateSpecies::initialise_arrays()
{
	nm_trans_ind.resize(sp.nm_trans.size());
	Li_nm_trans.resize(sp.nm_trans.size());
	Li_ie.resize(sp.ind_eff_group.size());
	
	auto G = sp.tra_gl.size();
	Li_markov_pop.resize(G);
		
	trans_num_f.resize(G);
	tnum_mean_st_f.resize(G);
	for(auto tr = 0u; tr < G; tr++){
		trans_num_f[tr].resize(T);
		tnum_mean_st_f[tr].resize(T);
	}
	
	pop_data_num.resize(sp.pop_data.size());
	pop_trans_data_num.resize(sp.pop_trans_data.size());
	
	obs_eqn_value.resize(sp.obs_eqn.size());
}


/// Initialises species ready for simulation be be performed
void StateSpecies::simulate_init(const InitCondValue &initc_val)
{
	auto G = sp.tra_gl.size();
	trans_num.clear(); trans_num.resize(G);
	tnum_mean_st.clear(); tnum_mean_st.resize(G);
	
	cpop_st.clear();
	
	ie_sampler_init();  // Samples individual effects
	
	markov_eqn_vari.clear();
	markov_eqn_vari.resize(N);
	for(auto e = 0u; e < N; e++){
		auto &me_vari = markov_eqn_vari[e];
	
		me_vari.time_vari = sp.markov_eqn[e].time_vari;
		
		if(me_vari.time_vari == true){
			me_vari.dt = details.dt; 
			me_vari.div.resize(T);
		}		
		else{
			me_vari.dt = UNSET;
			me_vari.div.resize(1);
		}			
		
		me_vari.indfac_sum = UNSET;
		me_vari.value = UNSET;
	}
	
	end_nm_event.clear();
	trig_div.clear();
	trig_div.resize(T);
		
	type = sp.type;

	cpop_init = initc_val.cpop;
	
	individual.clear();
	
	/// Adds any individuals from the data
	for(const auto &ind : sp.individual){
		if(ind.ev.size() == 0) emsg("No events");
		
		auto c_init = UNSET;
		
		auto i = individual.size();
		for(const auto &ev : ind.ev){
			auto t = ev.t;
			if(t == details.t_start){
				if(ev.type != ENTER_EV) emsg("SHould be enter");
				c_init = ev.c;
			}
			else{
				TrigEvent trig; trig.type = ev.type; trig.i = i; trig.t = t;
				
				switch(ev.type){
				case ENTER_EV:
					trig.c = ev.c; trig.trg = UNSET;
					c_init = UNSET;
					break;
					
				case LEAVE_EV:
					trig.c = UNSET; trig.trg = UNSET;
					break;
				
				case MOVE_EV:
					trig.c = ev.c; trig.trg = ev.cl;
					break;
					
				default: emsg("event type not recog"); break;
				}	
				insert_trigger_event(trig);
			}
		}
		add_individual(ind.name,c_init,false);
	}
	
	ie_Amatrix_sampler_init();
	sample_ie_Amatrix();	
	
	switch(type){
	case POPULATION:
		cpop = cpop_init;	
		break;

	case INDIVIDUAL:
		{
			cpop.clear();
			cpop.resize(sp.comp_gl.size(),0);
		
			string pre = ""; if(sp.individual.size() > 0) pre = "UO ";
			auto n = 0;
			for(auto c = 0u; c < sp.comp_gl.size(); c++){
				if(cpop_init[c] > 0 && sp.fix_effect.size() > 0){
					string st = "";
					for(const auto &fe : sp.fix_effect){
						if(st != "") st += ","; 
						st += fe.name;
					}
					emsg("Cannot add unspecified individuals for a model with fixed effect(s) '"+st+"'"); 
				}
				
				for(auto i = 0u; i < cpop_init[c]; i++){
					auto name = get_new_ind_name(pre,n+1);
					add_individual(name,c,false);
					n++;
				}
			}
		}
		break;
	}
	
	for(auto f = 0u; f < sp.fix_effect.size(); f++) set_exp_fe(f);
}


/// Gets a new unique name for an individual
string StateSpecies::get_new_ind_name(string pre, unsigned int num) const
{
	auto loop = 0u;
	do{
		string name;
		if(loop == 0) name = pre+"Ind. "+tstr(num); 
		else name = pre+"Ind. "+tstr(num)+"("+alphabet[loop-1]+")";
		
		auto i = 0u; while(i < individual.size() && individual[i].name != name) i++;
		if(i == individual.size()) return name;
		loop++;
	}while(loop < 27);
	
	emsg("Could not get individual name");
	return "";
}


/// Adds an individual to the state
void StateSpecies::add_individual(string name, unsigned int c, bool from_source)
{
	Individual ind;
	ind.name = name;
	ind.cinit = c;
	ind.c = c;
	ind.cl_ev.resize(sp.ncla);
	
	ind.ie = sample_ie();
	set_exp_ie(ind);
	
	if(from_source == true){
		if(sp.fix_effect.size() > 0) emsg("Cannot have fixed affects with a source");
	}
	else{
		auto i = individual.size();
		for(auto f = 0u; f < sp.fix_effect.size(); f++){
			const auto &fe = sp.fix_effect[f];
			const auto &X = fe.X_vector;
			
			if(i >= X.ind_list.size()) emsg("Individual '"+name+"' does not have a design matrix X specified for fixed effect '"+fe.name+"'");
			
			if(X.ind_list[i] != name) emsg("Names do not agree");
			ind.X.push_back(X.value[i]);
			ind.exp_fe.push_back(UNSET);
		}
	}
	
	individual.push_back(ind);
	
	if(individual.size() > details.individual_max){
		emsg("Maximum number of "+tstr(details.individual_max)+" individuals reached");
	}
}


/// Sets up the exponential of the fixed effects
vector <double> StateSpecies::set_exp_fe(unsigned int f)
{	
	const auto &fe = sp.fix_effect[f];
		
	const auto &par = param[fe.th];

	double pval;
	if(par.variety == CONST_PARAM) pval = par.value[0].value;
	else{
		auto th2 = par.param_vec_ref[0];
		if(th2 == UNSET) emsg("Should not be unset1");		
		pval = param_val[th2];
	}
	
	vector <double> store;

	auto sum = 0.0;
	for(auto &ind : individual){
		store.push_back(ind.exp_fe[f]);
		
		auto X = ind.X[f];
		if(X != UNSET){
			auto val = exp(X*pval);
			ind.exp_fe[f] = val;
			sum += val;
		}
	}
	
	auto fac = individual.size()/sum;

	for(auto &ind : individual){
		auto X = ind.X[f];
		if(X == UNSET) ind.exp_fe[f] = 1;
		else ind.exp_fe[f] *= fac;
	}
	
	return store;
}


/// Restores the values for exp_fe
void StateSpecies::set_exp_fe_restore(unsigned int f, const vector <double> &store)
{
	auto j = 0u;
	for(auto &ind : individual){
		ind.exp_fe[f] = store[j]; j++;
	}
}


/// Sets the values for exponentiated individual effects
vector <double> StateSpecies::set_exp_ie(Individual &ind) const
{
	const auto &indeff = sp.ind_effect;
	const auto &iegs = ind_eff_group_sampler;
	
	auto N = indeff.size();
	
	if(ind.ie.size() != N) emsg("Sizes do not agree");
	
	vector <double> store;
	
	ind.exp_ie.resize(N);
	for(auto i = 0u; i < indeff.size(); i++){
		const auto &ie = indeff[i];
		auto var = iegs[ie.index].omega[ie.num][ie.num];
		store.push_back(ind.exp_ie[i]);
		ind.exp_ie[i] = exp(ind.ie[i]-0.5*var);
	}
	
	return store;
}


/// Recalcualted exp_ie for a given individual effect
vector <double> StateSpecies::recalculate_exp_ie(unsigned int ie)
{
	auto &ind_eff = sp.ind_effect[ie];
	auto g = ind_eff.index;
	auto j = ind_eff.num;
	const auto &iegs = ind_eff_group_sampler[g];
	auto var = iegs.omega[j][j];
	
	vector <double> store;
	for(auto &ind : individual){
		store.push_back(ind.exp_ie[ie]);
		ind.exp_ie[ie] = exp(ind.ie[ie]-0.5*var);
	}
	
	return store;
}


/// Recalcualted exp_ie given omega has changed
void StateSpecies::recalculate_exp_ie_restore(unsigned int ie, const vector <double> &store)
{
	auto j = 0u;
	for(auto &ind : individual){
		ind.exp_ie[ie] = store[j]; j++;
	}
}


/// Initialises sampler used to sample individual effects
void StateSpecies::ie_sampler_init()
{
	ind_eff_group_sampler.clear();
	ind_eff_group_sampler.resize(sp.ind_eff_group.size());
	
	for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
		const auto &ieg = sp.ind_eff_group[i];
		auto &iegs = ind_eff_group_sampler[i];
		auto N = ieg.list.size();
		iegs.omega.resize(N); iegs.omega_Z.resize(N); iegs.omega_inv.resize(N);
		for(auto j = 0u; j < N; j++){
			iegs.omega[j].resize(N); iegs.omega_Z[j].resize(N); iegs.omega_inv[j].resize(N);
		}
		
		calculate_omega(i);
	}
}



/// 
void StateSpecies::ie_Amatrix_sampler_init()
{
	for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
		const auto &ieg = sp.ind_eff_group[i];
		auto &iegs = ind_eff_group_sampler[i];
		
		const auto &Aieg = ieg.A_matrix;
		if(Aieg.set == true){
			auto I = individual.size();
			vector <unsigned int> map(I);
			for(auto i = 0u; i < I; i++){
				auto name = individual[i].name;
				auto ii = find_in(Aieg.ind_list,name);
				if(ii == UNSET) emsg("Cannot find individual '"+name+"' is the Amatrix");
				map[i] = ii;
			}
			
			auto &A = iegs.A;
			A.resize(I);
			for(auto j = 0u; j < I; j++){
				A[j].resize(I);
				for(auto i = 0u; i < I; i++){
					A[j][i] = Aieg.value[map[j]][map[i]];
				}
			}
			iegs.A_Z = calculate_cholesky(A);
			auto A_inv = invert_matrix(A);
			iegs.A_inv = A_inv;
			iegs.A_inv_diag.resize(I);
			iegs.A_inv_nonzero.resize(I);
			for(auto j = 0u; j < I; j++){
				iegs.A_inv_diag[j] = A_inv[j][j];
				for(auto i = 0u; i < I; i++){
					if(A_inv[j][i] != 0){
						NZElement nze; nze.i = i; nze.value = A_inv[j][i];
						iegs.A_inv_nonzero[j].push_back(nze);
					}
				}
			}				
		}
	}
}
			

/// Calculates the omega matrix for individual effect
vector <double> StateSpecies::calculate_omega(unsigned int g)
{
	auto &iegs = ind_eff_group_sampler[g];
	const auto &ieg = sp.ind_eff_group[g];
	auto N = ieg.list.size();
	
	vector <double> store;
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			store.push_back(iegs.omega[j][i]);
			store.push_back(iegs.omega_Z[j][i]);
			store.push_back(iegs.omega_inv[j][i]);
		}
	}
	
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			auto &par = param[ieg.omega[j][i]];
			if(par.variety == CONST_PARAM) iegs.omega[j][i] = par.value[0].value;
			else{
				auto th2 = par.param_vec_ref[0]; if(th2 == UNSET) emsg("Should not be unset2");
				iegs.omega[j][i] = param_val[th2];
			}
		}
	}
	
	//print("Omega",iegs.omega);
	
	/// Converts correlations to variances
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			if(j != i) iegs.omega[j][i] *= sqrt(iegs.omega[i][i]*iegs.omega[j][j]);
		}
	}
	
	iegs.omega_Z = calculate_cholesky(iegs.omega);
	iegs.omega_inv = invert_matrix(iegs.omega);
	
	return store;
}


/// Restore the values in omega
void StateSpecies::calculate_omega_restore(unsigned int g, const vector <double> &store) 
{
	auto &iegs = ind_eff_group_sampler[g];
	const auto &ieg = sp.ind_eff_group[g];
	auto N = ieg.list.size();
	
	auto k = 0u;
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			iegs.omega[j][i] = store[k]; k++;
			iegs.omega_Z[j][i] = store[k]; k++;
			iegs.omega_inv[j][i] = store[k]; k++;
		}
	}
}


/// Samples the individual effects 
vector <double> StateSpecies::sample_ie() const
{
	auto N = sp.ind_effect.size(); // The number of individual effects
	
	vector <double> vec(N,UNSET);
	
	auto flag = false;
	for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
		const auto &ieg = sp.ind_eff_group[i];
		
		const auto &iegs = ind_eff_group_sampler[i];
		
		if(ieg.A_matrix.set == false){
			auto samp = sample_mvn(iegs.omega_Z);
		
			if(samp.size() != ieg.list.size()) emsg("Wrong size2");
			
			for(auto j = 0u; j < ieg.list.size(); j++){
				auto k = ieg.list[j].index;
				if(k >= N) emsg("Out of range5");
				if(vec[k] != UNSET) emsg("Should not be set");
				vec[k] = samp[j];
			}
		}
		else flag = true;
	}
	
	if(flag == false){
		for(auto i = 0u; i < N; i++){
			if(vec[i] == UNSET) emsg("SHould not be unset3");
		}
	}
	
	return vec;
}


/// Samples the individual effects with A matrix
void StateSpecies::sample_ie_Amatrix()
{
	auto flag = false;
	for(auto i = 0u; i < sp.ind_eff_group.size(); i++){
		const auto &ieg = sp.ind_eff_group[i];
		
		const auto &iegs = ind_eff_group_sampler[i];
		
		if(ieg.A_matrix.set == true){
			flag = true;
			
			auto N = ieg.list.size();
			auto I = individual.size();
			auto M = N*I;
			
			vector <double> vec(M);
			for(auto i = 0u; i < M; i++) vec[i] = normal_sample(0,1);
			
			for(auto j = 0u; j < M; j++){
				auto j_A = (unsigned int)(j/N);
				auto j_O = j%N;
				
				auto sum = 0.0;
				for(auto i = 0u; i < M; i++){
					auto i_A = (unsigned int)(i/N);
					auto i_O = i%N;
				
					sum += iegs.A_Z[j_A][i_A]*iegs.omega_Z[j_O][i_O]*vec[i];
				}
				
				individual[j_A].ie[ieg.list[j_O].index] = sum;
			}				
		}
	}
	
	for(auto &ind : individual){		
		set_exp_ie(ind);
		for(auto i = 0u; i < sp.ind_effect.size(); i++){
			if(ind.ie[i] == UNSET) emsg("Should be set1");
		}

		if(flag == true) set_exp_ie(ind);
	}
}


/// Returns the individual-based factor which acts on Markov equation
double StateSpecies::get_indfac(const Individual &ind, const MarkovEqn &mar_eqn) const
{
	auto val = 1.0;
	for(auto ie : mar_eqn.ind_eff_mult) val *= ind.exp_ie[ie];
	for(auto fe : mar_eqn.fix_eff_mult) val *= ind.exp_fe[fe];
	
	if(false){
		cout << ind.name << "   ";
		cout << eqn[mar_eqn.eqn_ref].te << " ";
		for(auto ie : mar_eqn.ind_eff_mult) cout << "mult IE " << ie << endl;
		for(auto fe : mar_eqn.fix_eff_mult) cout << "mult FE " << fe << endl;
		cout << val << "val\n";
	}
	
	return val;
}



/// Updates population-based species
void StateSpecies::update_population_based(unsigned int ti, bool stoc, const vector <double> &popnum)
{
	auto t = timepoint[ti];
	auto ddt = timepoint[ti+1] - t;
		
	auto tnum_mean = calculate_tnum_mean(ti,popnum,cpop,param_val,ddt);

	auto tnum = sample_trans_num(tnum_mean,stoc);
	//cout << ti << " " << tnum[0] << " " << tnum_mean[0] << " transnum\n";
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		tnum_mean_st[tr].push_back(tnum_mean[tr]);
		trans_num[tr].push_back(tnum[tr]);
	}
	
	cpop_st.push_back(cpop);
	
	update_cpop(cpop,tnum);
}


/// Updates the compartmental populations based on transition numbers
void StateSpecies::update_cpop(vector <double> &cpop, const vector <double> &tnum) const
{
	for(auto i = 0u; i < sp.tra_gl.size(); i++){
		const auto &tra = sp.tra_gl[i];
		auto num = tnum[i];
	
		if(tra.i != SOURCE) cpop[tra.i] -= num;
		if(tra.f != SINK) cpop[tra.f] += num;
	}
}


/// Samples transition number given transition rates
vector <double> StateSpecies::sample_trans_num(const vector <double> &tnum_mean, bool stochastic) const
{
	auto N = tnum_mean.size();
	vector <double> tnum(N);
	
	for(auto tr = 0u; tr < N; tr++){
		if(stochastic == true) tnum[tr] = poisson_sample(tnum_mean[tr]);	
		else tnum[tr] = tnum_mean[tr];
	}
	
	return tnum;
}


/// Calculates the rate for different transitions
vector <double> StateSpecies::calculate_tnum_mean(unsigned int ti, const vector <double> &popnum, const vector <double> &cpop, const vector <double> &param_val, double dt) const
{
	auto N = sp.tra_gl.size();
	vector <double> tnum_mean(N);
	
	for(auto i = 0u; i < N; i++){
		tnum_mean[i] = calculate_tnum_mean(ti,i,popnum,cpop,param_val,dt);
	}
	
	return tnum_mean;
}

/// Calculates the rate for different transitions
double StateSpecies::calculate_tnum_mean(unsigned int ti, unsigned int i, const vector <double> &popnum, const vector <double> &cpop, const vector <double> &param_val, double dt) const
{
	const auto &tr = sp.tra_gl[i];
	if(tr.type != EXP_RATE) emsg("Should be exponential");
	
	const auto &eq = eqn[tr.dist_param[0].eq_ref];
	auto rate = eq.calculate(ti,popnum,param_val,spline_val);
	
	if(rate < 0){
		for(auto val:param_val) cout << val << " val\n";
		
		emsg("The transition rate for '"+tr.name+"' has become negative."+check_prior(eq));
	}
		
	auto ci = tr.i;
	double pop;
	if(ci == SOURCE){ pop = 1.0/eq.source_tr_gl.size();}  
	else{ pop = cpop[ci]; if(pop < 0) pop = 0;}
	
	return pop*rate*dt;
}


/// Generates an A matrix (for testing)
void StateSpecies::generate_A()
{
	ofstream fout("Testing/M.csv");
	for(auto i = 0u; i < individual.size(); i++){
		if(i != 0) fout << ",";
		fout << "\"" << individual[i].name << "\""; 
	}		
	fout << endl;
	
	for(auto j = 0u; j < individual.size(); j++){
		for(auto i = 0u; i < individual.size(); i++){
			if(i != 0) fout << ",";
			if(i == j) fout << 1;
			else{
				if((unsigned int)(i/10) == (unsigned int)(j/10)) fout << 0.5;
				else fout << 0;
			}
		}
		fout << endl;
	}
}


/// Performs a MBP
void StateSpecies::mbp(double sim_prob, vector < vector <double> > &popnum_t)
{
	cpop_st_f.clear();
	
	auto cpop = cpop_st[0];
	//auto popnum = popnum_t[0];
	
	auto N = sp.tra_gl.size();
				
	vector <int> num(N);
	for(auto ti = 0u; ti < T; ti++){
		const auto &popnum = popnum_t[ti];
	
		cpop_st_f.push_back(cpop);
	
		auto ddt = dtimepoint[ti];
	
		for(auto tr = 0u; tr < N; tr++){
			auto tnum_mean_i = tnum_mean_st[tr][ti];
			auto tnum_mean_f = calculate_tnum_mean(ti,tr,popnum,cpop,param_val,ddt);
		
			auto num_i = trans_num[tr][ti];

			unsigned int num_f;	
			if(sim_prob == 0 || ran() > sim_prob){     // Does a MBP
				if(tnum_mean_f > tnum_mean_i){
					num_f = num_i + poisson_sample(tnum_mean_f-tnum_mean_i);
				}
				else{
					if(tnum_mean_i == tnum_mean_f) num_f = num_i;
					else{
						num_f = binomial_sample(tnum_mean_f/tnum_mean_i,num_i);
					}
				}
			}
			else{                                      // Does a simulation
				num_f = poisson_sample(tnum_mean_f);	
			}
			
			num[tr] = num_f;
			tnum_mean_st_f[tr][ti] = tnum_mean_f; 
			trans_num_f[tr][ti] = num_f;
		}
		
		/// Updates the populations
		if(ti < T-1){
			auto &popnum_new = popnum_t[ti+1];
			
			for(auto po : pop_affect) popnum_new[po] = popnum[po];
		
			for(auto tr = 0u; tr < N; tr++){
				auto n = num[tr];
				if(n != 0){
					const auto &tra = sp.tra_gl[tr];
					auto ci = tra.i;
					if(ci != SOURCE){
						cpop[ci] -= n;
						for(auto po : sp.comp_gl[ci].pop_ref_simp) popnum_new[po] -= n;
					}
					
					auto cf = tra.f;
					if(cf != SINK){
						cpop[cf] += n;
						for(auto po : sp.comp_gl[cf].pop_ref_simp) popnum_new[po] += n;
					}
				}
			}
		}
	}
	
	trans_num = trans_num_f;
	tnum_mean_st = tnum_mean_st_f;
	cpop_st = cpop_st_f;
}


/// Performs a MBP
void StateSpecies::mbp_accept(double &like_ch, const vector < vector <double> > &popnum_t)
{
	like_ch -= sum(Li_markov_pop);
	for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
		auto &Li = Li_markov_pop[tr];
		auto &tn = trans_num[tr];
		auto &tnm = tnum_mean_st[tr];
		for(auto ti = 0u; ti < T; ti++){
			tnm[ti] = calculate_tnum_mean(ti,tr,popnum_t[ti],cpop_st[ti],param_val,dtimepoint[ti]);
			Li[ti] = poisson_probability(tn[ti],tnm[ti]);
		}
	}
	like_ch += sum(Li_markov_pop);
}


