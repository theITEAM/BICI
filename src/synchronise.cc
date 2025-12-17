// Synchronise chains

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "synchronise.hh"
#include "mpi.hh"

/// Synchronises proposals across chains
void synchronise_proposal(unsigned int s, vector <Chain> &chain, Mpi &mpi)
{
	if(false) cout << mpi.core;
	
	if(s != UNSET){
		auto nburnin = chain[0].get_nburnin();
		if(s > nburnin) return;
		if(s != nburnin){
			auto step = nburnin/4;
			if(step == 0) step = 1;
			if(s%step != step/2) return;
		}
	}
	
	// Combines together the correlation matrices from all chains
	
	vector < vector <double> > all_samp;
	
	auto nchain = chain.size()*num_core();
	auto n = chain[0].cor_matrix.n;

	auto step = nchain;
	if(step == 0) step = 1;
	auto N = n/step;
		
	for(auto ch = 0u; ch < chain.size(); ch++){
		const auto &cha = chain[ch];
		const auto &cor_m = cha.cor_matrix;
		
		for(auto i = step-1; i < n; i += step){
			all_samp.push_back(cor_m.samp[i]);
		}
		
		//auto M = cor_m.calculate_cor_matrix();
		//print_matrix("mat before",M);
	}

#ifdef USE_MPI 
	mpi.transfer_samp(all_samp);
#endif
		
	vector < vector <double> > samp_new;

	if(op()){
		// Reorders in burning in sequence
		for(auto i = 0u; i < N; i++){
			for(auto ch = 0u; ch < nchain; ch++){
				samp_new.push_back(all_samp[ch*N+i]);
			}
		}
	}
	
#ifdef USE_MPI 
	mpi.distribute(samp_new);
#endif
	
	for(auto &cha : chain){
		auto &cor_m = cha.cor_matrix;
		cor_m.init();
	
		for(const auto &sa : samp_new){
			cor_m.add_sample2(sa,cor_m.n*FRAC_COR_UNUSED);
		}
		
		//auto M = cor_m.calculate_cor_matrix();
		//print_matrix("mat after",M);
		
		cha.join_proposal_update();
	}

	// Synchronise proposals
	
	//print_prop("BEFORE PROP",chain);
	
	vector <PropInfo> prop_info;
	
	for(const auto &cha : chain){
		auto ps_list = cha.get_prop_info();
		for(const auto &ps : ps_list) prop_info.push_back(ps);
	}
	
#ifdef USE_MPI 
	mpi.transfer_prop_info(prop_info);
#endif

	vector <PropInfo> prop_info_av;
	
	if(op()){ // Averages quantities
		vector <bool> done(prop_info.size(),false);
		
		for(auto i = 0u; i < prop_info.size(); i++){
			if(done[i] == false){
				auto &pi = prop_info[i];
				auto av = pi.value;
				auto vec_av = pi.vec;
				
				vector <unsigned int> list;
				for(auto j = i+1; j < prop_info.size(); j++){
					if(done[j] == false){
						const auto &pi2 = prop_info[j];
						if(pi2.type == pi.type && equal_vec(pi2.id,pi.id)){					
							av += pi2.value;
							for(auto k = 0u; k < vec_av.size(); k++) vec_av[k] += pi2.vec[k];
							done[j] = true;
							list.push_back(j);
						}
					}
				}
				
				auto ntot = list.size()+1;
				if(ntot != nchain){
					emsg("prop info nchain");
				}
				
				av /= ntot;
				for(auto k = 0u; k < vec_av.size(); k++) vec_av[k] /= ntot;

				pi.value = av;
				pi.vec = vec_av;
				prop_info_av.push_back(pi);
			}
		}
	}
	
#ifdef USE_MPI 
	mpi.distribute(prop_info_av);
#endif

	for(auto &cha : chain){
		cha.set_prop_info(prop_info_av);
	}

	//print_prop("AFTER PROP",chain);
}

				
/// Prints quantities used for proposals
void print_prop(string te, const vector <Chain> &chain)
{
	cout << te << endl;
	for(auto ch = 0u; ch < chain.size(); ch++){
		const auto &cha = chain[ch];
		for(auto i = 0u; i < cha.proposal.size(); i++){
			const auto &pro = cha.proposal[i];
			if(pro.on){
				switch(pro.type){	
				case PARAM_PROP:
					cout << pro.name << " " << pro.si << endl;
					break;
				default: break;
				}
			}
		}
	}		
}


/// Prints information from a prop_info
string print_prop_info(const vector <PropInfo> &prop_info, const Model &model)
{
	stringstream ss;
	for(const auto &pi : prop_info){
		ss << prop_type_str(pi.type) << ": ";
		switch(pi.type){
		case PARAM_PROP:
			for(auto th : pi.id){
				ss << model.param_vec[th].name << ",";
			}
			break;
			
		default: break;
		}
		
		ss << endl;
	}
	
	return ss.str();
}


