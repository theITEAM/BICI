/// Implements local population-based MCMC proposals

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

#include "proposal.hh"
#include "ind_ev_sampler.hh"
#include "utils.hh"
#include "matrix.hh"

/// Initialises 2D sampler
void Sampler2D::init(unsigned int _N, unsigned int _N2)
{
	N = _N; N2 = _N2;
	
	unsigned int range = N2/40;
	if(range > 10) range = 10;
	
	for(int dti = -2*range; dti <= (int)(2*range); dti++){
		KernelAdd ka;
		ka.dti = dti;
		ka.w = exp(normal_probability(dti,0,range));
		kernel.push_back(ka);
	}

	num.resize(N);    
	num_sum.resize(N);
	marg_sum.resize(N);
	for(auto i = 0u; i < N; i++){
		num[i].resize(N2,100);  
		num_sum[i].resize(N2); 
	}
	
	setup();
}


/// Sets up 2D sampler
void Sampler2D::setup()
{
	auto msum = 0.0;
	for(auto i = 0u; i < N; i++){
		auto sum = 0.0;
		for(auto j = 0u; j < N2; j++){
			sum += num[i][j];
			num_sum[i][j] = sum;
		}
		msum += sum;
		marg_sum[i] = msum;
	}
	S = msum;
}


/// Samples from 2D sampler
Samp Sampler2D::sample(double &dprob) const
{
	auto ra = ran();
	auto z = ra*S;

	auto i = 0u; while(i < N && z > marg_sum[i]) i++;
	if(i == N) emsg("Cannot sample i1");
	
	if(i > 0) z -= marg_sum[i-1];
	
	const auto &ns = num_sum[i]; 
	auto j = 0u; while(j < N2 && z > ns[j]) j++;
	if(j == N2) emsg("Cannot sample j");
	
	Samp samp; samp.i = i; samp.j = j;
	dprob -= log(num[i][j]/S);
	
	return samp;
}


/// Prints the sampler
void Sampler2D::print() const
{
	for(auto i = 0u; i < N; i++){
		cout << i << ": ";
		for(auto j = 0u; j < N2; j++){
			cout << int(num[i][j]) << ",";
		}
		cout << " num" << endl<< endl;
	}
	
	for(auto i = 0u; i < N; i++){
		cout << i << ": ";
		for(auto j = 0u; j < N2; j++){
			cout << int(num[i][j]) << ",";
		}
		cout << " num sum" << endl << endl;
	}
	
	for(auto i = 0u; i < N; i++){
		cout << i << ": " << marg_sum[i] << " marg" << endl;
	}
	
	cout << S << " S" << endl;
}


/// Updates the 2D sampler
void Sampler2D::update(unsigned int i, unsigned int j, int sign)
{
	num[i][j] += sign;
	
	auto &nsum = num_sum[i];
	for(auto jj = j; jj < N2; jj++) nsum[jj] += sign;
	
	for(auto ii = i; ii < N; ii++) marg_sum[ii] += sign;
	S += sign;
	if(S < 0) emsg("GG");
}


/// Check the sampler is correctly updated
void Sampler2D::check()
{
	auto num_sum_store = num_sum;
	auto marg_sum_store = marg_sum;
	auto S_store = S;
	
	setup();
	
	if(dif(num_sum_store,num_sum,DIF_THRESH)){
		emsg("problem with num sum");
	}
	
	if(dif(marg_sum_store,marg_sum,DIF_THRESH)){
		emsg("problem with marg sum");
	}
	
	if(dif(S_store,S,DIF_THRESH)){
		emsg("S");
	}
}


/// Initialises sample for compartments
void Sampler::init(unsigned int _N)
{
	N = _N;
	num.resize(N,5); 
	num_sum.resize(N);
	setup();
}


/// Sets up the sampler
void Sampler::setup()
{
	auto sum = 0.0;
	for(auto i = 0u; i < N; i++){
		sum += num[i];
		num_sum[i] = sum;
	}
	S = sum;
}


/// Samples from the sampler
unsigned int Sampler::sample(double &dprob) const
{
	auto z = ran()*S;
	auto i = 0u; while(i < N && z > num_sum[i]) i++;
	if(i == N) emsg("Cannot sample i2 ");

	dprob -= log(num[i]/S);
	
	return i;
}


/// Updates transition sampler
void Sampler::update(unsigned int i, int sign)
{
	num[i] += sign;
	
	for(auto ii = i; ii < N; ii++) num_sum[ii] += sign;
	S += sign;
}


/// Checks the sampler is correctly updated
void Sampler::check()
{
	auto num_sum_store = num_sum;
	auto S_store = S;
	
	setup();

	if(dif(num_sum_store,num_sum,DIF_THRESH)){
		emsg("problem with num sum");
	}
	
	if(dif(S_store,S,DIF_THRESH)){
		emsg("S");
	}
}


/// Returns the local window size
unsigned int Proposal::get_win() const
{
	return (unsigned int)(loc_samp.win);
}


/// Updates the window such that it is close to the optimum 
void Proposal::update_win(unsigned int num, unsigned int num_opt)
{
	if(!burn_info.on) return;
	
	if(num > 10*num_opt) num = 10*num_opt;
	if(num < 0.1*num_opt) num = 0.1*num_opt;
	auto f = WIN_UPDATE_ADAPT;
	auto win = f*loc_samp.win + (1-f)*loc_samp.win*(num_opt/(num+1));
	
	if(win < 1) win = 1;
	auto T = model.ntimepoint;
	if(win > 1+T/10) win = 1+T/10;
	loc_samp.win = win;
}


/// Adds or removes local changes to population-based model
void Proposal::pop_add_rem_local(State &state)
{
	const auto &tr_list = loc_samp.tr_list;
	auto T = state.T;
	auto L = tr_list.size();
	auto &ssp = state.species[p_prop];
	auto &tn = ssp.trans_num;
	
	auto win = get_win();
	
	// Randomly selects a window into which events are added or removed
	auto tmin = (unsigned int)((T-win+1)*ran());

	vector <Sampler> samp(L);	
	for(auto j = 0u; j < L; j++){
		samp[j].init(win);
		for(auto ti = 0u; ti < win; ti++) samp[j].num[ti] = tn[loc_samp.tr_list[j]][tmin+ti];
		samp[j].setup();
	}
		
	auto min = LARGE; for(auto j = 0u; j < L; j++){ if(samp[j].S < min) min = samp[j].S;}
	update_win(min,ADDREM_NUM_OPTIMUM);
	
	for(auto loop = 0u; loop < LOOP_ADD_REM_LOCAL; loop++){
		if(ran() < 0.5){                                // Add events
			vector <LocalChange> local_change;
		
			auto dprob = 0.0;
			for(auto j = 0u; j < L; j++){
				auto ti = (unsigned int)(win*ran());
				dprob += log(((samp[j].num[ti]+1)/(samp[j].S+1))/(1.0/win));
				LocalChange lc; 
				lc.type = ADD_EVENT; lc.tr = tr_list[j]; lc.ti = tmin+ti;
				local_change.push_back(lc);
			}
			
			auto like_ch = state.calculate_local_change(p_prop,local_change,1);
		
			auto al = calculate_al(like_ch,dprob);
			
			ntr++;
			if(ran() < al){
				nac++;	
				
				for(auto j = 0u; j < L; j++) samp[j].update(local_change[j].ti-tmin,1);
				state.add_like(like_ch);
			}
			else{
				state.calculate_local_change(p_prop,local_change,-1);
			}
		}
		else{                                           // Remove events
			auto j = 0u; while(j < L && samp[j].S > 0) j++;
			if(j == L){
				auto dprob = 0.0;
				
				vector <LocalChange> local_change;
				for(auto j = 0u; j < L; j++){
					auto ti = samp[j].sample(dprob);
					LocalChange lc; 
					lc.type = REMOVE_EVENT; lc.tr = tr_list[j]; lc.ti = tmin+ti;
					local_change.push_back(lc);
				}
				dprob += L*log(1.0/win);
				
		
				auto like_ch = state.calculate_local_change(p_prop,local_change,1);
			
				auto al = calculate_al(like_ch,dprob);
			
				ntr++;
				if(ran() < al){
					nac++;	
					for(auto j = 0u; j < L; j++) samp[j].update(local_change[j].ti-tmin,-1);
					state.add_like(like_ch);
				}
				else{ 
					state.calculate_local_change(p_prop,local_change,-1);
				}
			}
		}		
	}
	
	if(false){
		state.check("hh");
		for(auto j = 0u; j < L; j++) samp[j].check();
	}
}

	
/// Moves local events in population-based model
void Proposal::pop_move_local(State &state)
{
	auto T = state.T;
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	auto &tn = ssp.trans_num;
	auto N = sp.tra_gl.size();
	
	auto win = (unsigned int)(LOCAL_MOVE_RANGE*T);
	
	auto tmin = (unsigned int)((T-win+1)*ran());

	Sampler2D samp;
	samp.init(N,win);
	for(auto tr = 0u; tr < N; tr++){
		for(auto ti = 0u; ti < win; ti++) samp.num[tr][ti] = tn[tr][tmin+ti];
	}
	samp.setup();
	
	if(samp.S == 0) return;
			
	auto temp = 0.0;
	for(auto loop = 0u; loop < LOOP_MOVE_LOCAL; loop++){
		auto po = samp.sample(temp);
		auto tr = po.i, ti = tmin+po.j;
		
		auto ti_new = tmin + (unsigned int)(ran()*win);
		
		vector <LocalChange> local_change;
		LocalChange lc; 
		lc.type = REMOVE_EVENT; lc.tr = tr; lc.ti = ti;
		local_change.push_back(lc);
		
		lc.type = ADD_EVENT; lc.tr = tr; lc.ti = ti_new;
		local_change.push_back(lc);
		
		auto like_ch = state.calculate_local_change(p_prop,local_change,1);
		
		auto al = calculate_al(like_ch,0);
		
		ntr++;
		if(ran() < al){
			nac++;
			samp.update(tr,ti-tmin,-1);
			samp.update(tr,ti_new-tmin,1);
			state.add_like(like_ch);
		}
		else{ 
			state.calculate_local_change(p_prop,local_change,-1);
		}
	}
	
	if(false){
		state.check("hh");
		for(auto tr = 0u; tr < N; tr++){
			for(auto ti = 0u; ti < win; ti++){
				if(dif(samp.num[tr][ti],tn[tr][tmin+ti],DIF_THRESH)) emsg("dif prob");
			}
		}
	}
}

	
/// Makes changes to init cond in population-based model
void Proposal::pop_ic_local(State &state)
{
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	auto &tn = ssp.trans_num;
	auto &cnum = ssp.init_cond_val.cnum;
	auto N = sp.tra_gl.size();
		
	auto win = get_win();
	
	auto &tr_ic = loc_samp.tr_ic;
	
	if(burn_info.on){                                 // Updates the sampler
		for(auto tr = 0u; tr < N; tr++){
			for(auto ti = 0u; ti < win; ti++) tr_ic.num[tr] += tn[tr][ti];
		}
		tr_ic.setup();
	}
	
	Sampler2D samp;
	samp.init(N,win);
	for(auto tr = 0u; tr < N; tr++){
		for(auto ti = 0u; ti < win; ti++) samp.num[tr][ti] = tn[tr][ti];
	}
	samp.setup();
		
	update_win(samp.S,IC_NUM_OPTIMUM);
	
	for(auto loop = 0u; loop < LOOP_IC_LOCAL; loop++){
		if(ran() < 0.5){                                // Add event and change initial condition	
			auto dprob = 0.0;
			auto tr = tr_ic.sample(dprob);
			
			auto cf = sp.tra_gl[tr].f; 
			if(cf == UNSET || cnum[cf] > 0){
				auto ti = (unsigned int)(ran()*win);
				dprob += log(((samp.num[tr][ti]+1)/(samp.S+1))/(1.0/win));
				
				vector <LocalChange> local_change;
				LocalChange lc; 
				lc.type = ADD_IC; lc.tr = tr; lc.ti = -1;
				local_change.push_back(lc);
				
				lc.type = ADD_EVENT; lc.tr = tr; lc.ti = ti;
				local_change.push_back(lc);
				
				auto like_ch = state.calculate_local_change(p_prop,local_change,1);
				
				auto al = calculate_al(like_ch,dprob);
			
				ntr++;
				if(ran() < al){ 
					nac++;
					samp.update(tr,ti,1);
					state.add_like(like_ch);
				}
				else{ 
					state.calculate_local_change(p_prop,local_change,-1);
				}
			}
		}
		else{                                           // Remove event and change initial condition
			if(samp.S > 0){
				auto dprob = 0.0;
				auto sa = samp.sample(dprob);
				auto tr = sa.i, ti = sa.j;
				
				auto ci = sp.tra_gl[tr].i; 
				if(ci == UNSET || cnum[ci] > 0){
					dprob += log((tr_ic.num[tr]/tr_ic.S)*(1.0/win));
				
					vector <LocalChange> local_change;
					LocalChange lc; 
					lc.type = REMOVE_IC; lc.tr = tr; lc.ti = -1;
					local_change.push_back(lc);
				
					lc.type = REMOVE_EVENT; lc.tr = tr; lc.ti = ti;
					local_change.push_back(lc);
				
					auto like_ch = state.calculate_local_change(p_prop,local_change,1);
				
					auto al = calculate_al(like_ch,dprob);
			
					ntr++;
					if(ran() < al){ 
						nac++;
						samp.update(tr,ti,-1);
						state.add_like(like_ch);
					}
					else{ 
						state.calculate_local_change(p_prop,local_change,-1);
					}
				}
			}
		}
	}
	
	if(false){
		for(auto tr = 0u; tr < N; tr++){
			for(auto ti = 0u; ti < win; ti++){
				if(dif(samp.num[tr][ti],tn[tr][ti],DIF_THRESH)) emsg("Wrong");
			}
		}
		state.check("hh");
		samp.check();
		tr_ic.check();
	}
}


/// Makes changes to end in population-based model
void Proposal::pop_end_local(State &state)
{
	auto T = state.T;
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	auto &tn = ssp.trans_num;
	auto N = sp.tra_gl.size();
		
	auto win = get_win();
	
	auto &tr_ic = loc_samp.tr_ic;
	
	auto ti_shift = T-win;
	
	if(burn_info.on){                                 // Updates the sampler
		for(auto tr = 0u; tr < N; tr++){
			for(auto ti = 0u; ti < win; ti++) tr_ic.num[tr] += tn[tr][ti+ti_shift];
		}
		tr_ic.setup();
	}
	
	Sampler2D samp;
	samp.init(N,win);
	for(auto tr = 0u; tr < N; tr++){
		for(auto ti = 0u; ti < win; ti++) samp.num[tr][ti] = tn[tr][ti+ti_shift];
	}
	samp.setup();
		
	update_win(samp.S,IC_NUM_OPTIMUM);
	
	for(auto loop = 0u; loop < LOOP_END_LOCAL; loop++){
		if(ran() < 0.5){                                // Add event 
			auto dprob = 0.0;
			auto tr = tr_ic.sample(dprob);
			
			auto ti = ti_shift+(unsigned int)(ran()*win);	
			dprob += log(((samp.num[tr][ti-ti_shift]+1)/(samp.S+1))/(1.0/win));
				
			vector <LocalChange> local_change;
			LocalChange lc; 	
			lc.type = ADD_EVENT; lc.tr = tr; lc.ti = ti;
			local_change.push_back(lc);
				
			auto like_ch = state.calculate_local_change(p_prop,local_change,1);
			
			auto al = calculate_al(like_ch,dprob);
		
			ntr++;
			if(ran() < al){ 
				nac++;
				samp.update(tr,ti-ti_shift,1);
				state.add_like(like_ch);
			}
			else{ 
				state.calculate_local_change(p_prop,local_change,-1);
			}
		}
		else{                                           // Remove event 
			if(samp.S > 0){
				auto dprob = 0.0;
				auto sa = samp.sample(dprob);
				auto tr = sa.i, ti = ti_shift+sa.j;
				
				dprob += log((tr_ic.num[tr]/tr_ic.S)*(1.0/win));
			
				vector <LocalChange> local_change;
				LocalChange lc; 
				lc.type = REMOVE_EVENT; lc.tr = tr; lc.ti = ti;
				local_change.push_back(lc);
			
				auto like_ch = state.calculate_local_change(p_prop,local_change,1);
			
				auto al = calculate_al(like_ch,dprob);
		
				ntr++;
				if(ran() < al){ 
					nac++;
					samp.update(tr,ti-ti_shift,-1);
					state.add_like(like_ch);
				}
				else{ 
					state.calculate_local_change(p_prop,local_change,-1);
				}
			}
		}
	}
	
	if(false){
		for(auto tr = 0u; tr < N; tr++){
			for(auto ti = 0u; ti < win; ti++){
				if(dif(samp.num[tr][ti],tn[tr][ti+ti_shift],DIF_THRESH)) emsg("Wrong");
			}
		}
		state.check("hh");
		samp.check();
		tr_ic.check();
	}
}


/// Adds or removes single events (with potentially making corresponding change to IC)
void Proposal::pop_single_local(State &state)
{
	auto pl = false;//true;
	
	const auto T = model.ntimepoint-1;
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	auto &tn = ssp.trans_num;
	auto &cnum = ssp.init_cond_val.cnum;
	auto N = sp.tra_gl.size();
	
	// This is a sampler for adding new events
	auto &tr_samp = loc_samp.tr_samp;
		
	if(burn_info.on){ // Updates the sampler
		vector <unsigned int> tr_num(N);
		for(auto tr = 0u; tr < N; tr++){	
			auto &num = tr_samp.num[tr];

			for(auto ti = 0u; ti < T; ti++){
				auto n = tn[tr][ti];
				if(n > 0){
					for(const auto &ka : tr_samp.kernel){
						int tii = ti + ka.dti;
						if(tii >= 0 && tii < int(T)) num[tii] += n*ka.w;
					}
				}
			}
		}
		tr_samp.setup();
	}
	
	// Setup removal sampler
	Sampler2D tr_sel;
	tr_sel.init(N,T);
	tr_sel.num = tn;
	tr_sel.setup();

	for(auto loop = 0u; loop < LOOP_IC_LOCAL; loop++){
		if(sp.init_cond.type == INIT_POP_DIST && ran() < 0.5){  // Changes IC and add/rem event
			if(ran() < 0.5){                               // Add event with initical condition change
				auto dprob = 0.0;
				auto samp = tr_samp.sample(dprob);
				auto tr = samp.i;
				auto ti = samp.j;
				
				auto cf = sp.tra_gl[tr].f; 
				if(cf == UNSET || cnum[cf] > 0){
					dprob += log((tr_sel.num[tr][ti]+1)/(tr_sel.S+1));
						
					vector <LocalChange> local_change;
					LocalChange lc; 
					lc.type = ADD_IC; lc.tr = tr; lc.ti = -1;
					local_change.push_back(lc);
				
					lc.type = ADD_EVENT; lc.tr = tr; lc.ti = ti;
					local_change.push_back(lc);
						
					auto like_ch = state.calculate_local_change(p_prop,local_change,1);
					
					auto al = calculate_al(like_ch,dprob);
					
					ntr++;
					if(ran() < al){ 
						nac++;	
						tr_sel.update(tr,ti,1);
						state.add_like(like_ch);
					}
					else{ 
						state.calculate_local_change(p_prop,local_change,-1);
					}
				}
			}
			else{                                         // Remove event and change initial condition
				if(tr_sel.S > 0){
					auto dprob = 0.0;
					auto samp = tr_sel.sample(dprob);
				
					auto tr = samp.i;	
					auto ti = samp.j;
				
					auto ci = sp.tra_gl[tr].i; 
					if(ci == UNSET || cnum[ci] > 0){
						dprob += log(tr_samp.num[tr][ti]/tr_samp.S);
						
						vector <LocalChange> local_change;
						LocalChange lc; 
						lc.type = REMOVE_IC; lc.tr = tr; lc.ti = -1;
						local_change.push_back(lc);
				
						lc.type = REMOVE_EVENT; lc.tr = tr; lc.ti = ti;
						local_change.push_back(lc);
					
						auto like_ch = state.calculate_local_change(p_prop,local_change,1);
				
						auto al = calculate_al(like_ch,dprob);
				
						ntr++;
						if(ran() < al){ 
							nac++;
							tr_sel.update(tr,ti,-1);
							state.add_like(like_ch);
						}
						else{ 
							state.calculate_local_change(p_prop,local_change,-1);
						}
					}
				}
			}
		}
		else{                                           // Just add/rem event
			if(ran() < 0.5){ // Add event 
				if(pl) cout << "add" << endl;
				
				auto dprob = 0.0;
				auto samp = tr_samp.sample(dprob);
				auto tr = samp.i;
				auto ti = samp.j;
				
				dprob += log((tr_sel.num[tr][ti]+1)/(tr_sel.S+1));
					
				vector <LocalChange> local_change;
				LocalChange lc; 	
				lc.type = ADD_EVENT; lc.tr = tr; lc.ti = ti;
				local_change.push_back(lc);
					
				auto like_ch = state.calculate_local_change(p_prop,local_change,1);
					
				auto al = calculate_al(like_ch,dprob);
				
				ntr++;
				if(ran() < al){ 
					nac++;	
					tr_sel.update(tr,ti,1);
					state.add_like(like_ch);
				}
				else{ 
					state.calculate_local_change(p_prop,local_change,-1);
				}
			}
			else{                                         // Remove event and change initial condition
				if(tr_sel.S > 0){
					if(pl) cout << "remove" << endl;
					
					auto dprob = 0.0;
					auto samp = tr_sel.sample(dprob);
				
					auto tr = samp.i;	
					auto ti = samp.j;
					
					dprob += log((tr_samp.num[tr][ti])/(tr_samp.S));
					
					vector <LocalChange> local_change;
					LocalChange lc; 
					lc.type = REMOVE_EVENT; lc.tr = tr; lc.ti = ti;
					local_change.push_back(lc);
				
					auto like_ch = state.calculate_local_change(p_prop,local_change,1);
			
					auto al = calculate_al(like_ch,dprob);
			
					ntr++;
					if(ran() < al){ 
						nac++;
						tr_sel.update(tr,ti,-1);
						state.add_like(like_ch);
					}
					else{ 
						state.calculate_local_change(p_prop,local_change,-1);
					}
				}
			}
		}
			
		if(pl) state.check("Proposal local");
	}
	
	if(false){
		if(dif(tr_sel.num,tn,DIF_THRESH)) emsg("Not agree");
	
		state.check("hh");
		tr_sel.check();
	}
}


/// Makes changes to init cond in population-based model
void Proposal::pop_ic(State &state)
{
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	auto &cnum = ssp.init_cond_val.cnum;
	auto N = sp.comp_gl.size();
	auto &comp_ic = loc_samp.comp_ic;
	
	if(burn_info.on){ // Updates the sampler
		for(auto c = 0u; c < N; c++) comp_ic.num[c] += cnum[c];
		comp_ic.setup();
	}
	
	// Setup removal sampler
	Sampler samp;
	samp.init(N);
	for(auto c = 0u; c < N; c++) samp.num[c] = cnum[c];
	samp.setup();
		
	for(auto loop = 0u; loop < LOOP_IC; loop++){
		if(ran() < 0.5){                                // Add event with initical condition change
			auto dprob = 0.0;
			auto c = comp_ic.sample(dprob);
			if(!sp.comp_gl[c].erlang_hidden){
				dprob += log((samp.num[c]+1)/(samp.S+1));
						
				vector <LocalChange> local_change;
				LocalChange lc; 
				lc.type = ADD_C_IC; lc.tr = c; lc.ti = -1;
				local_change.push_back(lc);
					
				auto like_ch = state.calculate_local_change(p_prop,local_change,1);
					
				auto al = calculate_al(like_ch,dprob);
			
				ntr++;
				if(ran() < al){
					nac++;	
					samp.update(c,1);
					state.add_like(like_ch);
				}
				else{ 
					state.calculate_local_change(p_prop,local_change,-1);
				}
			}
		}
		else{                                           // Remove event and change initial condition
			if(samp.S > 0){
				auto dprob = 0.0;
				auto c = samp.sample(dprob);
				if(!sp.comp_gl[c].erlang_hidden){
					dprob += log(comp_ic.num[c]/comp_ic.S);
						
					vector <LocalChange> local_change;
					LocalChange lc; 
					lc.type = REMOVE_C_IC; lc.tr = c; lc.ti = -1;
					local_change.push_back(lc);
			
					auto like_ch = state.calculate_local_change(p_prop,local_change,1);
			
					auto al = calculate_al(like_ch,dprob);
			
					ntr++;
					if(ran() < al){ 
						nac++;
						samp.update(c,-1);
						state.add_like(like_ch);
					}
					else{ 
						state.calculate_local_change(p_prop,local_change,-1);
					}
				}
			}
		}
	}
	
	if(false){
		for(auto c = 0u; c < N; c++){
			if(dif(samp.num[c],double(cnum[c]),DIF_THRESH)) emsg("Not agree");
		}

		state.check("hh");
		samp.check();
		comp_ic.check();
	}
}


/// Swaps c in init cond for population-based model
void Proposal::pop_ic_swap(State &state)
{
	const auto &sp = model.species[p_prop];
	auto &ssp = state.species[p_prop];
	auto &cnum = ssp.init_cond_val.cnum;
	auto N = sp.comp_gl.size();
	auto &comp_ic = loc_samp.comp_ic;
	
	if(burn_info.on){ // Updates the sampler
		for(auto c = 0u; c < N; c++) comp_ic.num[c] += cnum[c];
		comp_ic.setup();
	}
	
	// Setup removal sampler
	Sampler samp;
	samp.init(N);
	for(auto c = 0u; c < N; c++) samp.num[c] = cnum[c];
	samp.setup();
	
	if(samp.S == 0) return;
	
	for(auto loop = 0u; loop < LOOP_IC_SWAP; loop++){
		auto dprob = 0.0;
			
		auto ci = samp.sample(dprob);
		auto cf = comp_ic.sample(dprob);
		dprob += log(((samp.num[cf]+1)/(samp.S+1))*(comp_ic.num[ci]/comp_ic.S));
								
		vector <LocalChange> local_change;
		LocalChange lc; 
		lc.type = REMOVE_C_IC; lc.tr = ci; lc.ti = -1;
		local_change.push_back(lc);
		
		lc.type = ADD_C_IC; lc.tr = cf; lc.ti = -1;
		local_change.push_back(lc);
				
		auto like_ch = state.calculate_local_change(p_prop,local_change,1);
				
		auto al = calculate_al(like_ch,dprob);
		
		ntr++;
		if(ran() < al){
			nac++;	
			samp.update(ci,-1);
			samp.update(cf,1);
			state.add_like(like_ch);
		}
		else{ 
			state.calculate_local_change(p_prop,local_change,-1);
		}
	}
	
	if(false){
		for(auto c = 0u; c < N; c++){
			if(dif(samp.num[c],double(cnum[c]),DIF_THRESH)) emsg("Not agree");
		}

		state.check("hh");
		samp.check();
		comp_ic.check();
	}
}
