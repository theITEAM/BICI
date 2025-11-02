// This file contains output functions

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <algorithm> 
#include <iomanip>
 
using namespace std;

#include "output.hh"
#include "utils.hh"
//#include "state.hh"

/// Initialises the output 
Output::Output(const Model &model, const Input &input, Mpi &mpi) : model(model), mpi(mpi)
{
	timer.resize(OUTTIMER_MAX); for(auto &ti : timer) ti = 0;
	
	diagdir = "";	sampledir = ""; sampledir_rel = "";

	if(input.datadir != ""){
		switch(model.mode){
		case SIM: sampledir_rel = "output-sim"; break;
		case INF: sampledir_rel = "output-inf"; break;
		case PPC: sampledir_rel = "output-post-sim"; break;
		case MODE_UNSET: emsg("Option not recognised"); break;
		}
		
		sampledir = input.datadir+"/"+sampledir_rel;
		ensure_directory(sampledir);  
		
		if(model.mode == INF){
			diagdir = sampledir+"/diagnostic";
			ensure_directory(diagdir);                  // Creates the diagnostics directory
		}
		
		//if(op()) start_progess_file(diagdir+"/progress.txt");
	}
	
	if(op()){
		lines_raw = input.lines_raw;
		
		auto i = 0u;
		while(i < lines_raw.size())
		{
			auto st = trim(lines_raw[i]);
			
			auto spl = split(st,' ');
			auto command = spl[0];
			
			auto flag = false;
			
			if(begin_str(st,"# PROCESSED USING")) flag = true;

			if(command == "param-sim" || command == "state-sim" || command == "warning-sim" || st == "# OUTPUT SIMULATION"){
				if(model.mode == SIM) flag = true;
			}
			
			if(command == "param-inf" || command == "param-stats-inf" || command == "state-inf" || 
			   command == "diagnostics-inf" || command == "generation-inf" || 
				 command == "trans-diag-inf" || command == "warning-inf" || st == "# OUTPUT INFERENCE"){
				if(model.mode == INF) flag = true;
			}
			
			if(command == "param-post-sim" || command == "state-post-sim" || command == "warning-post-sim" || st == "# OUTPUT POSTERIOR SIMULATION"){
				if(model.mode == INF || model.mode == PPC) flag = true;
			}
	
			if(flag == true){
				if(st.length() >= 3 && st.substr(st.length()-3,3) == "\"[["){
					auto j = i+1;
					while(j < lines_raw.size() && !(lines_raw[j].length() >= 3 && lines_raw[j].substr(0,3) == "]]\"")) j++;
					if(j == lines_raw.size()) emsg("Problem removing output");
					lines_raw.erase(lines_raw.begin()+i,lines_raw.begin()+j+1); 
				}
				else{
					lines_raw.erase(lines_raw.begin()+i); 
				}
			}
			else i++;
		}
		
		while(lines_raw.size() > 0 && lines_raw[lines_raw.size()-1] == ""){
			lines_raw.erase(lines_raw.begin()+lines_raw.size()-1); 
		}
		
		// Removes more than double space
		auto j = 0u;
		while(j+2 < lines_raw.size()){
			if(trim(lines_raw[j]) == "" && trim(lines_raw[j+1]) == "" && trim(lines_raw[j+2]) == ""){
				lines_raw.erase(lines_raw.begin()+j);
			}
			else j++;
		}
		
		lines_raw.push_back("");
		switch(model.mode){
		case SIM: lines_raw.push_back("# OUTPUT SIMULATION"); break;
		case INF: lines_raw.push_back("# OUTPUT INFERENCE"); break;
		case PPC: lines_raw.push_back("# OUTPUT POSTERIOR SIMULATION"); break;
		case MODE_UNSET: break;
		}
		lines_raw.push_back("");
	}
};


/// Create a directory if it doesn't already exist
void Output::ensure_directory(const string &path) const 
{
	if(!op()) return;
	
	struct stat st;
	if (stat(path.c_str(), &st) == -1){  	            // Directory not found

#ifdef WINDOWS
		int ret = mkdir(path.c_str());
#else	
		int ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
		
		if(ret == -1){
			run_error("Error creating the directory '"+path+"'");
		}
	}
}


/// Checks that the file has been successfully open
void Output::check_open(ofstream &fout, string file) const 
{
	if(!fout) run_error("Could not write to file '"+file+"'");
}


/// Outputs a summary of the model and data to a file
void Output::summary(const Model &model) const
{
	if(diagdir == "") return;
	auto file = diagdir+"/model.txt";
	ofstream fout(file);
	check_open(fout,file);

	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		fout << "SPECIES: " << sp.name << "   type: ";
		switch(sp.type){
		case POPULATION: fout << "population-based"; break;
		case INDIVIDUAL: fout << "individual-based"; break;
		}

		fout << endl << endl;
		
		fout << "N: " << sp.N << endl;

		for(auto ci = 0u; ci < sp.ncla; ci++){
			const auto &cl = sp.cla[ci];
	
			fout << "Classification: " << cl.name << "  index: " << cl.index << "  coord: ";
			switch(cl.coord){
			case CARTESIAN: fout << "cartesian"; break;
			case LATLNG: fout << "latlng"; break;
			}
			fout << endl << endl;
			
			fout << "Compartments:" << endl;
			for(const auto &co : cl.comp){
				fout << co.name << "  branch: ";
				if(co.branch == true) fout << "true"; else fout << "false";
				
				if(sp.trans_tree == true && sp.infection_cl == ci){
					fout <<   "  infected: ";
					if(co.infected == COMP_INFECTED) fout << "true"; else fout << "false";
				}
				
				switch(cl.coord){
				case CARTESIAN: fout << "  x: " << co.x << "  y: " << co.y; break;
				case LATLNG: fout << "  lat: " << co.lat << "  lng: " << co.lng; break;
				}
				fout << endl;
			}		

			fout << endl;
	
			fout << "Transitions:" << endl;
			for(const auto &tr : cl.tra){
				fout << tr.name << "   ";
				if(tr.i == UNSET) fout << "+"; else fout << cl.comp[tr.i].name;

				fout << " -> ";
				if(tr.f == UNSET) fout << "-"; else fout << cl.comp[tr.f].name;

				fout << "    " << transtype_text(tr.type);
				
				fout << "    ";
				switch(tr.variety){
				case NORMAL: fout << "normal"; break;
				case SOURCE_TRANS: fout << "source"; break;
				case SINK_TRANS: fout << "sink"; break; 
				}
				
				fout << "   bp_set: "; if(tr.bp_set == true) fout << "true"; else fout << "false";
				fout << "   branch: "; if(tr.branch == true) fout << "true"; else fout << "false";
			
				fout << endl;
			}
			fout << endl;
		}
	
		if(false){
			fout << "Global Compartments:" << endl;
			for(const auto &co_gl : sp.comp_gl){
				fout << co_gl.name;
				if(sp.trans_tree == true){
					fout << "  infected: ";
					if(co_gl.infected == true) fout << "true"; else fout << "false";
				}
				fout << endl;
				for(auto cl = 0u; cl < sp.ncla; cl++){
					const auto &cla = sp.cla[cl];
					
					const auto &tlg = co_gl.tra_leave_group[cl];
					
					fout << cla.name << " ";
					if(tlg.branch == true) fout << "branch"; else fout << "no branch";
					fout << ": ";
					for(auto tgl : tlg.tr_list){
						const auto &tra = sp.tra_gl[tgl];
						fout << tra.name;
						
						if(tlg.branch == true){
							fout << " BP=";
							if(tra.bp_set == BP_DERIVED) fout << "* ";
							else fout << model.eqn[tra.bp.eq_ref].te_raw << " ";
						}
						fout << ", ";
					}	
					fout << endl;
				}
				
				fout << "pop affect: ";
				for(auto po : co_gl.pop_ref_simp) fout << model.pop[po].name << ", ";
				fout << endl;
			}

			fout << "Global Transitions:" << endl;
			for(const auto &tr_gl : sp.tra_gl){
				fout << tr_gl.name << "     ";
				if(tr_gl.i == UNSET) fout << "+";
				else fout << sp.comp_gl[tr_gl.i].name;
				fout << " -> ";
				if(tr_gl.f == UNSET) fout << "-";
				else fout << sp.comp_gl[tr_gl.f].name;
		
				fout << "  " << transtype_text(tr_gl.type) << "("; 
				for(auto j = 0u; j <  tr_gl.dist_param.size(); j++){
					if(j != 0) fout << ",";
					fout << model.eqn[tr_gl.dist_param[j].eq_ref].te_raw;
				}
				fout << ")";
		
				if(sp.trans_tree == true){
					switch(tr_gl.infection.type){
					case TRANS_INFECTION: fout << "  infection"; break;
					case TRANS_RECOVERY: fout << "  recovery"; break;
					case TRANS_INF_UNCHANGE: break;
					}
				}
		
				if(tr_gl.branch == true){
					if(tr_gl.bp_set == BP_DERIVED) fout << "     BP = *";
					else fout << "    BP = " << model.eqn[tr_gl.bp.eq_ref].te_raw;
				}
			
				fout << "   time vari: "; if(tr_gl.time_vari == true) fout << "true"; else fout << "false";
			
				fout << endl;
			}
			
			fout << "NM trans:" << endl;
			for(const auto &nmt : sp.nm_trans){
				fout << transtype_text(nmt.type);
				fout << " all_branches:" << nmt.all_branches;
				fout << " precalc_nm_rate:" << nmt.precalc_nm_rate;
				if(nmt.bp_eq != UNSET){
					fout << " bp=";
					if(nmt.bp_eq == BP_FROM_OTHERS) fout << "*"; else fout << model.eqn[nmt.bp_eq].te;
					fout << " bp_other=";
					for(auto e : nmt.bp_other_eq) fout << model.eqn[e].te << ",";
					fout << " bp_all=";
					for(auto e : nmt.bp_all_eq) fout << model.eqn[e].te << ",";
				}
				fout << " dist=";
				for(auto e : nmt.dist_param_eq_ref) fout << model.eqn[e].te << ",";
				fout << endl;
			}
			
			fout << "NM trans incomp:" << endl;
			for(const auto &nmti : sp.nm_trans_incomp){
				for(auto m : nmti.nmtrans_ref) fout << m << ",";
				fout << endl;
			}
		}
		
		fout << endl;

		for(auto ci = 0u; ci < sp.ncla; ci++){
			const auto &cl = sp.cla[ci];
			
			fout << "Local Swap for " << cl.name << ":" << endl; 
			for(const auto &sr : cl.swap_rep){
				fout << sr.name << endl;
			}
			fout << endl;
		}
	
		fout << "DATA:" << endl;
		
		for(const auto &ds : sp.source){
			switch(ds.cname){
			case INIT_POP_SIM: fout << "Init. pop sim"; break;
			case ADD_POP_SIM: fout << "Add pop sim"; break;
			case REMOVE_POP_SIM: fout << "Remove pop sim"; break;
			case ADD_IND_SIM: fout << "Add ind sim"; break;
			case REMOVE_IND_SIM: fout << "Remove ind sim"; break;
			case MOVE_IND_SIM: fout << "Move ind sim"; break; 
			case ADD_POP_POST_SIM: fout << "Add pop post-sim"; break;
			case REMOVE_POP_POST_SIM: fout << "Remove pop post-sim"; break;
			case ADD_IND_POST_SIM: fout << "Add ind post-sim"; break;
			case REMOVE_IND_POST_SIM: fout << "Remove ind post-sim"; break;
			case MOVE_IND_POST_SIM: fout << "Move ind post-sim"; break; 
			case INIT_POP: fout << "Init pop"; break; 
			case ADD_POP: fout << "Add pop"; break; 
			case REMOVE_POP: fout << "Remove pop"; break;
			case ADD_IND: fout << "Add ind"; break; 
			case REMOVE_IND: fout << "Remove ind"; break;
			case MOVE_IND: fout << "Move ind"; break; 
			case COMP_DATA: fout << "Comp data"; break; 
			case TRANS_DATA: fout << "Trans data"; break;
			case TEST_DATA: fout << "Test data"; break; 
			case POP_DATA: fout << "Pop data"; break; 
			case POP_TRANS_DATA: fout << "Pop trans data"; break; 
			case GENETIC_DATA: fout << "Genetic data"; break;
			case PARAM_MULT: fout << "Parameter multiply"; break;
			default: fout << "Unknown data type"; break;
			}
	
			fout << "  "; 
			const auto &tab = ds.table;
			fout << tab.file << "  ncol: " << tab.ncol << "  nrow: " << tab.nrow;
			fout << endl;
		}
		fout << endl;
		
		fout << "NM_TRANS" << endl;
		for(const auto &nmt : sp.nm_trans){
			fout << transtype_text(nmt.type) << "  ";
		
			if(nmt.bp_eq != UNSET){
				if(nmt.bp_eq == BP_FROM_OTHERS) fout << "BP=*  ";
				else fout << "BP=" << model.eqn[nmt.bp_eq].te_raw << "  ";
			}
			for(auto eq : nmt.dist_param_eq_ref) fout << model.eqn[eq].te_raw << "  ";
			fout << endl;
		}
		fout << endl;
	
		fout << "FIXED EFFECTS" << endl;
		for(const auto &fe : sp.fix_effect){
			fout << fe.name << "   param:" << fe.th << endl;
		}
		fout << endl;
		
		fout << "INDIVIDUAL EFFECT GROUPS" << endl;
		for(const auto &ieg : sp.ind_eff_group){
			for(const auto &li : ieg.list) fout << li.name;
			fout << endl;
			fout << "Omega:" << endl;
			auto N = ieg.list.size();
			const auto &par = model.param[ieg.th];
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i < N; i++){
					fout << par.name << "_" << ieg.list[j].name << "," << ieg.list[i].name << "  ";
				}
				fout << endl;
			}
			fout << "Markov_eqn_ref:";
			for(auto e : ieg.markov_eqn_ref) fout << e << " ";
			fout << endl;
			fout << "nm_trans_ref:";
			for(auto e : ieg.nm_trans_ref) fout << e << " ";
			fout << endl;
		}
		fout << endl;
	}
	
	if(model.genetic_data.on){
		fout << "GENETIC DATA:" << endl;
		fout << "# SNP=" << model.genetic_data.nSNP << endl;
		fout << "# obs=" << model.genetic_data.obs.size() << endl;
		fout << endl;
	}
	
	fout << "POPULATIONS:" << endl;
	for(const auto &po : model.pop){
		const auto &sp = model.species[po.sp_p];
		
		fout << po.name << " in " << sp.name << ": ";
		for(auto k = 0u; k < po.ind_eff_mult.size(); k++){
			fout << "IE " << sp.ind_effect[po.ind_eff_mult[k]].name << "  "; 
		}
		
		for(auto k = 0u; k < po.fix_eff_mult.size(); k++){
			fout << "FE " << sp.fix_effect[po.fix_eff_mult[k]].name << "  "; 
		}
		
		for(const auto &te : po.term){
			fout << sp.comp_gl[te.c].name << "," << te.w << "  ";
		}
		
		if(false){
			fout << " ME: ";
			for(const auto &mer : po.markov_eqn_ref){
				const auto &sp = model.species[mer.p];
				auto eq = sp.markov_eqn[mer.e].eqn_ref;
				fout << model.eqn[eq].te_raw << "  ";
			}
		}
		fout << endl;
	}
		
	fout << endl;
	
	fout << "INDIVIDUAL OBSERVATIONS" << endl;
	if(false){
		for(auto p = 0u; p < model.nspecies; p++){
			const auto &sp = model.species[p];
		
			fout << "SPECIES: " << sp.name << "   type: " << endl;
	
			for(const auto &ind : sp.individual){
				fout << ind.name << ":" << endl;
				for(const auto &ob : ind.obs){
					fout << "t=" << model.calc_t(ob.tdiv) << " ";
					switch(ob.type){
					case OBS_TRANS_EV: case OBS_SOURCE_EV: case OBS_SINK_EV:
						fout << "Transition: ";
						for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
							fout << sp.tra_gl[tr].name << "=";
							unsigned int oe;
							if(ob.time_vari) oe = ob.obs_eqn_ref[tr];
							else oe = sp.obs_trans[ob.ref].obs_eqn_ref[tr];
							fout << model.eqn[sp.obs_eqn[oe]].te_raw << ", ";
						}
						break;
						
					case OBS_COMP_EV:
						{
							fout << "Compartment: ";
							const auto &claa = sp.cla[ob.cl];
							for(auto c = 0u; c < claa.comp.size(); c++){
								fout << claa.comp[c].name << "=";
								fout << model.eqn[sp.obs_eqn[ob.obs_eqn_ref[c]]].te_raw << ", ";
							}
						}
						break;
						
					case OBS_TEST_EV:
						fout << ob.test_res << " ";
						fout << "Se=" << model.eqn[sp.obs_eqn[ob.Se_obs_eqn_ref]].te_raw << ", ";
						fout << "Sp=" << model.eqn[sp.obs_eqn[ob.Sp_obs_eqn_ref]].te_raw << ", ";
						break;
						
					default: emsg("Model output error 4"); break;
					}
					fout << endl;
				}
				fout << endl;
			}
		}	
	}
	else{
		fout << "Not output" << endl;
	}

	fout << "PARAMETERS" << endl;

	for(const auto &par : model.param){
		fout << par.name << "   ";
	
		switch(par.variety){
		case PRIOR_PARAM: fout << "prior"; break;
		case REPARAM_PARAM: fout << "reparam"; break;
		case DIST_PARAM: fout << "dist"; break;
		case CONST_PARAM: fout << "const"; break;
		case UNSET_PARAM: emsg("Error"); break;
		}
		
		fout << "  ";
		fout << "time_dep: "; if(par.time_dep == true) fout << "true"; else fout << "false";
		fout << endl;

		for(const auto &dep : par.dep){
			fout << "  index: " << dep.index << "  mult: " << dep.mult << " ";
			for(const auto &val : dep.list) fout << val << ",";	
			fout << endl;
		} 

		if(par.spline_info.on == true){
			fout << "  Spline: ";
			fout << "  ";
			fout << "smooth: "; if(par.spline_info.smooth == true) fout << "true"; else fout << "false";
			if(par.spline_info.smooth == true){
				switch(par.spline_info.smooth_type){
				case NORMAL_SMOOTH: fout << "Normal(" << par.spline_info.smooth_value << ")"; break;
				case LOG_NORMAL_SMOOTH: fout << "Log-Normal(" << par.spline_info.smooth_value << ")"; break;
				}
			}
			fout << "  times: ";
			for(auto tdiv : par.spline_info.knot_tdiv) fout << model.calc_t(tdiv) << ",";
			fout << endl;
		}

		switch(model.mode){
		case SIM:
			/*
			if(par.value.size() > 0){
				auto st = stringify(par.value,par.dep);
				fout << "  value: " << st;
			}
			
			
			fout << " parent: [";
			for(auto i = 0u; i < par.parent.size(); i++){
				if(i != 0) fout << ",";
				const auto &pp = par.parent[i];
				for(auto j = 0u; j < pp.size(); j++){
					if(j != 0) fout << "|";
					fout << model.param[pp[j].th].name << pp[j].index;
				}
			}			
			fout << "]  ";
			
			fout << " child: [";
			for(auto i = 0u; i < par.child.size(); i++){
				if(i != 0) fout << ",";
				const auto &pp = par.child[i];
				for(auto j = 0u; j < pp.size(); j++){
					if(j != 0) fout << "|";
					fout << model.param[pp[j].th].name << pp[j].index;
				}
			}			
			fout << "]  ";
			*/
			break;

		case INF: case PPC:
			break;

		default:  emsg("Model output error10"); break;
		}
		
		fout << endl;
		fout << endl;
	}
	
	if(false){
		fout << "PARAMETER VECTOR:" << endl;
		for(const auto &par : model.param_vec){
			fout << par.name;
		
			if(false){
			//if(true){
				fout << "   affect: " << endl;
				for(auto &al : par.affect_like){
					fout << print_affect_like(al);
				}
			}
			fout << endl;
		}
		fout << endl;
	}
		
	if(true){
		fout << "SPLINES:" << endl;
		for(const auto &spl : model.spline){
			fout << spl.name << endl;
		}
		fout << endl;
	}
	
	if(false){
		fout << "EQUATIONS:" << endl;
		for(const auto &eq : model.eqn){
			fout << eq.te_raw << "  ";
			if(eq.source_tr_gl.size() > 0){
				fout << "SOURCE TRANS: ";
				for(auto tr : eq.source_tr_gl) fout << model.species[eq.sp_p].tra_gl[tr].name <<", ";
			}
			fout << endl;
		}
		fout << endl;
	}
	
	if(false){
		fout << "DERIVED:" << endl;
		for(const auto &der : model.derive){
			fout << der.name << endl;
			for(const auto &eq : der.eq){
				fout << eq.te_raw << "   ref: " << eq.eq_ref << endl; 
			}
		}
		fout << endl;
	}
}


/// Outputs a summary of proprosals
void Output::prop_summary(string te) const
{
	if(!op()) return;
		
	if(debugging){
		auto file = "proposal.txt";
		ofstream fout(file);
		check_open(fout,file);
		fout << te;
		return;
	}

	if(diagdir == "") return;
	auto file = diagdir+"/proposal.txt";
	ofstream fout(file);
	check_open(fout,file);
	fout << te;
}


// Prints how the state properties are altered
string Output::print_affect_like(const AffectLike &al) const
{
	stringstream ss;
	
	switch(al.type){
	case OBS_EQN_AFFECT:
		{
			const auto &sp = model.species[al.num];
			ss << "AFFECT obs eqn for species " << sp.name;
		}
		break;
		
	case LIKE_UNOBS_TRANS_AFFECT:
		{
			const auto &sp = model.species[al.num];
			ss << "AFFECT like unobs trans for species " << sp.name << " equation " << model.eqn[sp.obs_trans_eqn[al.num2]].te;
		}
		break;
		
	case POP_DATA_CGL_TGL_AFFECT:
		{
			const auto &sp = model.species[al.num];
			ss << "AFFECT pop data cgl tgl for species " << sp.name;
		}
		break;
		
	case LIKE_OBS_IND_AFFECT:
		{
			ss << "AFFECT like obs ind for species " << model.species[al.num].name;
		}
		break;
		
	case LIKE_OBS_POP_AFFECT:
		{
			ss << "AFFECT like obs pop for species " << model.species[al.num].name;
		}
		break;
		
	case LIKE_OBS_POP_TRANS_AFFECT:
		{
			ss << "AFFECT like obs pop_trans for species " << model.species[al.num].name;
		}
		break;
		
	case MARKOV_POP_AFFECT:
		{
			ss << "AFFECT markov pop" << model.species[al.num].tra_gl[al.num2].name;
		}
		break;
		
	case LIKE_IE_AFFECT:
		{
			const auto &ieg = model.species[al.num].ind_eff_group[al.num2];
			ss << "AFFECT IE likelihood: ";
			for(auto li : ieg.list) ss << li.name << ", ";
		}
		break;
	
	case LIKE_INIT_COND_AFFECT:
		{
			ss << "AFFECT likelihood init cond: " <<  model.species[al.num].name;
		}
		break;
		
	case PRIOR_INIT_COND_AFFECT:
		{
			ss << "AFFECT prior init cond ";
		}
		break;
		
	case EXP_IE_AFFECT:
		{
			const auto &ind_eff = model.species[al.num].ind_effect[al.num2];
			ss << "AFFECT  exp ie: " <<  ind_eff.name;
		}
		break;
	
	case OMEGA_AFFECT:
		{
			const auto &ieg = model.species[al.num].ind_eff_group[al.num2];
			ss << "AFFECT omega: ";
			for(auto li : ieg.list) ss << li.name << ", ";
		}
		break;
				
	case NM_TRANS_AFFECT:
		{
			const auto &sp = model.species[al.num];
			const auto &nm = sp.nm_trans[al.num2];
			ss << "AFFECT Nonmarkovian trans affect: " << nm.name << " ";
			ss << model.eqn[nm.dist_param_eq_ref[0]].te_raw << endl;	
		}
		break;
		
	case NM_TRANS_BP_AFFECT:
		{
			const auto &nm = model.species[al.num].nm_trans[al.num2];
			ss << "AFFECT Nonmarkovian BP affect: ";
			if(nm.bp_eq == BP_FROM_OTHERS) ss << "*" << endl;
			else ss << model.eqn[nm.bp_eq].te_raw << endl;	
		}
		break;
		
	case NM_TRANS_INCOMP_AFFECT:
		{
			const auto &nm = model.species[al.num].nm_trans_incomp[al.num2];
			ss << "AFFECT Incomplete Nonmarkovian trans affect: " << nm.name << endl;	
		}
		break;
		
	case POP_AFFECT:
		{
			ss << "AFFECT Population affect" << endl;	
		}
		break;
		
	case SPLINE_PRIOR_AFFECT:
		{
			ss << "AFFECT Spline Prior affect " << model.spline[al.num].name << endl;	
		}
		break;
		
	case IEG_PRIOR_AFFECT:
		{
			const auto &iegr = model.ieg_ref[al.num];
			ss << "AFFECT IEG Prior affect " << model.species[iegr.p].ind_eff_group[iegr.i].name << endl;	
		}
		break;
		
	case PRIOR_AFFECT:
		{
			ss << "AFFECT Prior affect " << model.param_vec[al.num].name << endl;	
		}
		break;
	
	case DIST_AFFECT:
		{
			ss << "AFFECT Distribution affect " << model.param_vec[al.num].name << endl;	
		}
		break;
	
	case DIV_VALUE_AFFECT:
		{
			auto eq = model.species[al.num].markov_eqn[al.num2].eqn_ref;
			ss << "AFFECT Div Value " << model.eqn[eq].te_raw << endl;
		}
		break;
		
	case DIV_VALUE_NOPOP_AFFECT:
		{
			auto eq = model.species[al.num].markov_eqn[al.num2].eqn_ref;
			ss << "AFFECT Div Value nopop" << model.eqn[eq].te_raw << endl;
		}
		break;
		
	case DIV_VALUE_LINEAR_AFFECT:
		{
			ss << "AFFECT Div Value Linear ";
			{
				if(al.lin_form.factor_nopop_only) ss << "[factor/no-pop] ";
				else ss << "[full] ";
				string str = "";	
				for(const auto &lf : al.lin_form.list){
					str += model.eqn[lf.e].te_raw +", ";
				}
				ss << trunc(str,1000);
			}
			ss << endl;
		}
		break;
		
	case MARKOV_POP_NOPOP_AFFECT:
		{
			auto eq = model.species[al.num].tra_gl[al.num2].dist_param[0].eq_ref;
			ss << "AFFECT Markov pop nopop" << model.eqn[eq].te_raw << endl;
		}
		break;
		
	case MARKOV_POP_LINEAR_AFFECT:
		{
			ss << "AFFECT Markov pop linear ";
			{
				if(al.lin_form.factor_nopop_only) ss << "[factor/no-pop] ";
				else ss << "[full] ";
				string str = "";
				for(const auto &lf : al.lin_form.list){
					str += model.eqn[lf.e].te_raw +", ";
				}
				ss << trunc(str,1000);
			}
			ss << endl;
		}
		break;
		
	case INDFAC_INT_AFFECT:
		{
			auto &sp = model.species[al.num];
			ss << "AFFECT Indfac Ind  species" << sp.name << endl;
			
		}
		break;
		
	case EXP_FE_AFFECT:
		{
			const auto &fe = model.species[al.num].fix_effect[al.num2];
			ss << "AFFECT Set exp_fe for " << fe.name << endl;
		}
		break;
	
	case MARKOV_LIKE_AFFECT: 
		{
			auto eq = model.species[al.num].markov_eqn[al.num2].eqn_ref;
			ss << "AFFECT Markov Like " << model.eqn[eq].te_raw << endl; 
		}
		break;
		
	case LIKE_GENETIC_PROCESS_AFFECT:
		ss << "AFFECT Genetic Process Like " << endl;
		break;
		
		
	case GENETIC_VALUE_AFFECT:
		ss << "AFFECT Genetic Value " << endl;
		break;
		
	case LIKE_GENETIC_OBS_AFFECT:
		ss << "AFFECT Genetic Obs Like " << endl;
		break;
		
	case IIF_W_AFFECT:
		ss << "AFFECT IIF W" << endl;
		break;
		
	case POPNUM_IND_W_AFFECT:
		{
			const auto &pop = model.pop[al.num];
			ss << "AFFECT popnum ind w " << pop.name << endl;
		}
		break;
		
	case AFFECT_MAX: break;
	}
	
	if(al.list.size() != 0){
		switch(al.type){
		case POP_AFFECT:
			ss << "   Populations affected: ";	
			for(auto k : al.list) ss << model.pop[k].name << ","; 
			break;
		
		case OBS_EQN_AFFECT:
			ss << "   Obs Eqn affected: ";	
			for(auto i : al.list){
				const auto &eqn = model.eqn[model.species[al.num].obs_eqn[i]];
				ss << eqn.te_raw << ", ";
			}
			break;
			
		case LIKE_OBS_IND_AFFECT:
			ss << "   Ind affected: ";	
			if(al.list.size() ==  model.species[al.num].individual.size()) ss << "All";
			else{
				for(auto i : al.list) ss << model.species[al.num].individual[i].name << ",";
			}
			break;
				
		case LIKE_OBS_POP_AFFECT:
			ss << "   Pop data affected: ";	
			for(auto i : al.list) ss << i << ",";
			break;
							
		case LIKE_OBS_POP_TRANS_AFFECT:
			ss << "   Pop trans data affected: ";	
			for(auto i : al.list) ss << i << ",";
			break;
			
		default:	
			ss << "   Times affected: ";	
			
			if(al.list.size() == al.map.size()) ss << " All";
			else ss << model.str_time_range(al.list);	
			break;
		}
		ss << endl;	
	}

	return ss.str();	
}


/// Converts a type into text;
string Output::transtype_text(TransType type) const 
{
	switch(type){
	case EXP_RATE: return "Exp(rate)"; 
	case EXP_MEAN: return "Exp(mean)"; 
	case EXP_RATE_NM: return "Exp(rate) NM"; 
	case EXP_MEAN_NM: return "Exp(mean) NM"; 
	case GAMMA: return "Gamma";
	case ERLANG: return "Erlang"; 
	case LOG_NORMAL: return "Log-normal"; 
	case WEIBULL: return "Weibull";
	case PERIOD: return "Period";
	}
	return "Error";
}


/// Outputs a summary of the data to a file
void Output::data_summary(const Model &model) const
{
	if(!op() || diagdir == "") return;

	auto file = diagdir+"/data.txt";
	ofstream fout(file);
	check_open(fout,file);
	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		fout << "SPECIES: " << sp.name << endl;
		
		fout << "ENTER:" << endl;
		for(const auto &en : sp.enter){
			fout << en.name << " "<< model.calc_t(en.tdiv);
			if(en.c_set != UNSET){
				fout << "  " << sp.comp_gl[en.c_set].name << endl;
			}
			else{
				fout << endl;
				for(auto cl = 0u; cl < sp.ncla; cl++){
					const auto &claa = sp.cla[cl];
					const auto &en_cl = en.cla[cl];
					if(en_cl.c_set != UNSET){
						fout << " " << claa.comp[en_cl.c_set].name << endl;
					}
					else{
						fout << " ";
						auto flag = false;
						for(auto c = 0u; c < claa.ncomp; c++){
							auto te = en_cl.eqn[c].te;
							if(te != "0"){
								if(flag == true) fout << ","; 
								flag = true;
								fout << claa.comp[c].name << ":" << te;
							}
						}
					}
				}
			}
		}
	
		for(const auto &ind : sp.individual){
			fout << ind.name << ": ";
			for(const auto &ev : ind.ev){
				switch(ev.type){
				case ENTER_EV: fout << "enter [" << ind.enter_ref <<"]"; break;
				case LEAVE_EV: fout << "leave"; break;
				case MOVE_EV: fout << "move to " << sp.cla[ev.cl].comp[ev.move_c].name; break;
				case M_TRANS_EV: fout << "markov trans " << sp.cla[ev.cl].tra[ev.tr].name; break;	
				case NM_TRANS_EV: fout << "non-markov trans " << sp.cla[ev.cl].tra[ev.tr].name; 
				break;
				}
				fout << "," << model.calc_t(ev.tdiv) << "  ";
			}
			fout << endl;
			
			fout << print_obs_data(p,ind.obs);
		}
		fout << endl;
		
		fout << "Population filter:" << endl;
		for(const auto &pf : sp.pop_filter){
			fout << pf.name << endl;
			if(pf.time_vari == false){
				for(auto c : pf.c_nonzero){ 
					auto om = pf.comp_obs_mod_ref[c];
					const auto &eq = model.eqn[sp.obs_eqn[om]];
					fout << "  " << sp.comp_gl[c].name << " eq=" << eq.te_raw;
					fout << " obsmodref=" << om << endl;
				}
			}
		}
		
		fout << "Population data:" << endl;
		for(const auto &pd : sp.pop_data){
			fout << "t=" << model.calc_t(pd.tdiv) << "  value=" << pd.value;
			switch(pd.type){
			case NORMAL_OBS: fout << "  normal(sd=" << pd.obs_mod_val << ")"; break;
			case POISSON_OBS: fout << "  poisson()"; break;
			case NEGBIN_OBS: fout << "  neg-binomial(p=" << pd.obs_mod_val << ")"; break;
			}
		
			const auto &pf = sp.pop_filter[pd.ref];
			if(pf.time_vari == true){
				for(auto c : pf.c_nonzero){ 
					auto om = pd.comp_obs_mod_ref[c];
					const auto &eq = model.eqn[sp.obs_eqn[om]];
					fout << "  " << sp.comp_gl[c].name << " eq=" << eq.te_raw;
					fout << " obsmodref=" << om << endl;
				}
			}
		}
		fout << endl;
		
		fout << "Pop Trans filter:" << endl;
		for(const auto &ptf : sp.pop_trans_filter){
			fout << ptf.name << endl;
			if(ptf.time_vari == false){
				for(auto tr : ptf.tr_nonzero){ 
					auto om = ptf.trans_obs_mod_ref[tr];
					const auto &eq = model.eqn[sp.obs_eqn[om]];
					fout << "  " << sp.tra_gl[tr].name << " eq=" << eq.te_raw;
					fout << " obsmodref=" << om << endl;
				}
			}
		}
		
		fout << "Pop Trans data:" << endl;
		for(const auto &ptd : sp.pop_trans_data){
			fout << "start=" << model.calc_t(ptd.tdivmin);
			fout << "  end=" << model.calc_t(ptd.tdivmax);
			fout << "  value=" << ptd.value;
			
			switch(ptd.type){
			case NORMAL_OBS: fout << "  normal(sd=" << ptd.obs_mod_val << ")"; break;
			case POISSON_OBS: fout << "  poisson()"; break;
			case NEGBIN_OBS: fout << "  neg-binomial(p=" << ptd.obs_mod_val << ")"; break;
			}
			fout << endl;
			
			const auto &ptf = sp.pop_trans_filter[ptd.ref];
			if(ptf.time_vari == true){
				for(auto tr : ptf.tr_nonzero){ 
					auto om = ptd.trans_obs_mod_ref[tr];
					const auto &eq = model.eqn[sp.obs_eqn[om]];
					fout << "  " << sp.tra_gl[tr].name << " eq=" << eq.te_raw;
					fout << " obsmodref=" << om << endl;
				}
			}
		}
		fout << endl;
		
		fout << "Observation equation:" << endl;
		for(auto eq : sp.obs_eqn){
			fout << model.eqn[eq].te_raw << endl;
		}
	}
}


/// Prints a set of initial conditions
void Output::print_initc(const vector <InitCondValue> &initc) const
{
	if(com_op == true) return;
	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		const auto &ic = sp.init_cond;
	
		cout << "Initial conditions for '"+sp.name+"':" << endl;
		
		switch(ic.type){
		case INIT_POP_FIXED:
			cout << "Fixed:" << endl;
			for(auto c = 0u; c < sp.comp_gl.size(); c++){
				cout << sp.comp_gl[c].name << " " << ic.cnum[c] << endl;
			}
			break;
	
		case INIT_POP_DIST:
			if(ic.focal_cl != UNSET){
				auto fc = ic.focal_cl;
				cout << "Focal: " << sp.cla[fc].name << endl;
				for(auto c = 0u; c < sp.cla[fc].ncomp; c++){
					cout << "Prior for " << sp.cla[fc].comp[c].name << ": ";
					cout << get_prior_string(ic.comp_prior[c]) << endl;
				}
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(cl != fc){
						cout << "Alpha for population in " << sp.cla[cl].name << ": ";
						for(auto c = 0u; c < sp.cla[cl].ncomp; c++){
							const auto &co = sp.cla[cl].comp[c];
							if(!co.erlang_hidden) cout << ic.alpha_focal[cl][c] << ", "; 
						}
						cout << endl;
					}
				}
			}
			else{
				cout << "Total population prior: ";
				cout << get_prior_string(ic.pop_prior) << endl;
				cout << "Alpha values:" << endl;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					const auto &co = sp.comp_gl[c];
					if(!co.erlang_hidden) cout << sp.comp_gl[c].name << " "<< ic.alpha[c] << endl; 
				}			
			}
			break;
			
		case INIT_POP_NONE:
			cout << "No initial population set" << endl;
			return;
		}
	
		cout << "Value:" << endl;
		for(auto c = 0u; c < sp.comp_gl.size(); c++){
			cout << initc[p].cnum.size() << " " << c << " hhh" << endl;
			cout << sp.comp_gl[c].name << " " << initc[p].cnum[c] << endl;
		}
	}
	cout << endl;
}


/// Outputs a parameter sample 
string Output::trace_init() const 
{
	stringstream ss;
	
	ss << "State";
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		
		if(par.trace_output){
			auto pn = replace_arrow(par.name);
			pn = add_escape_char(pn);
			
			if(par.dep.size() == 0){		
				ss << ",\"" << pn << "\"";
			}
			else{
				if(model.is_symmetric(par)){  // Outputs covariance matrix
					auto L = par.dep[0].list.size();
					for(auto j = 0u; j < L; j++){  // Diagonal elements
						ss << ",\"" << pn << "_";
						ss << par.dep[0].list[j] << "," << par.dep[0].list[j];
						ss << "\"";
					}
					
					for(auto j = 0u; j < L; j++){  // Off-diagonal elements
						for(auto i = j+1; i < L; i++){	
							ss << ",\"" << model.exchange_omega(pn) << "_";
							ss << par.dep[0].list[j] << "," << par.dep[0].list[i];
							ss << "\"";
						}
					}
				}
				else{
					for(auto j = 0u; j < par.N; j++){
						ss << ",\"" << pn << "_";
						for(auto k = 0u; k < par.dep.size(); k++){
							const auto &dp = par.dep[k];
							auto m = (unsigned int)(j/dp.mult)%dp.list.size();
							if(k != 0) ss << ",";
							ss << dp.list[m];
						}
						ss << "\"";
					}
				}
			}
		}
	}
	
	for(auto i = 0u; i < model.derive.size(); i++){
		const auto &der = model.derive[i];
		if(der.time_dep == false){
			if(der.dep.size() == 0){
				ss << ",\"" << der.name << "\"";
			}
			else{
				for(auto j = 0u; j < der.eq.size(); j++){
					ss << ",\"" << der.name << "_";
					for(auto k = 0u; k < der.dep.size(); k++){
						const auto &dp = der.dep[k];
						auto m = (unsigned int)(j/dp.mult)%dp.list.size();
						if(k != 0) ss << ",";
						ss << dp.list[m];
					}
					ss << "\"";
				}
			}
		}
	}
	
	if(model.trans_tree){
		ss << ",N^origin,N^infected,N^mut-tree,N^mut-origin,N^unobs,t^root";
	}
	
	for(const auto &sp : model.species){
		if(sp.type == INDIVIDUAL) ss << ",N^" << sp.name << "-total";;
	}
	
	ss << ic_output_head();
	
	ss << ",L^markov,L^non-markov,L^ie,L^dist,L^obs,L^genetic-proc,L^genetic-obs,L^init";

	if(model.mode == INF) ss << ",Prior";
	
	ss << endl;
	
	return ss.str();
}


/// Outputs the parameter file columns for initial conditions
string Output::ic_output_head() const 
{
	stringstream ss;
	
	for(const auto &sp : model.species){
		const auto &ic = sp.init_cond;
		if(ic.type == INIT_POP_DIST){
			if(ic.focal_cl == UNSET){
				ss << ",N^" << sp.name;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					const auto &co = sp.comp_gl[c];
					if(co.erlang_hidden == false) ss << ",f^init_(" << co.name << ")";
				}
			}
			else{
				const auto &claa = sp.cla[ic.focal_cl];
				for(auto c = 0u; c < claa.comp.size(); c++){
					const auto co = claa.comp[c];
					if(!co.erlang_hidden) ss << ",N^init_(" << co.name << ")";
				}
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(cl != ic.focal_cl){
						const auto &claa = sp.cla[cl];
						for(auto c = 0u; c < claa.comp.size(); c++){
							const auto co = claa.comp[c];
							if(!co.erlang_hidden) ss << ",f^init_(" << co.name << ")";
						}
					}
				}
			}
		}
	}
	
	return ss.str();
}
	
	
/// Outputs the initial conditions
string Output::ic_output(const Particle &part) const
{
	stringstream ss;
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		const auto &ic = sp.init_cond;
		const auto &icv = part.species[p].init_cond_val;
		
		if(ic.type == INIT_POP_DIST){
			if(ic.focal_cl == UNSET){
				ss << "," << icv.N_total;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					const auto &co = sp.comp_gl[c];
					if(!co.erlang_hidden) ss << "," << icv.frac[c];
				}
			}
			else{
				const auto &claa = sp.cla[ic.focal_cl];
				for(auto c = 0u; c < claa.comp.size(); c++){
					const auto &co = claa.comp[c];
					if(!co.erlang_hidden) ss << "," << icv.N_focal[c];
				}
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(cl != ic.focal_cl){
						const auto &claa = sp.cla[cl];
						for(auto c = 0u; c < claa.comp.size(); c++){
							const auto &co = claa.comp[c];
							if(!co.erlang_hidden) ss << "," << icv.frac_focal[cl][c];
						}
					}
				}
			}
		}
	}

	return ss.str();
}


/// Outputs the burnin fraction (needed for anneal scan)	
void Output::set_output_burnin(double burnin_frac)
{
	if(!op()) return;
	
	for(auto j = 0u; j < lines_raw.size(); j++){
		auto st = lines_raw[j];
		auto spl = split(st,' ');
		if(spl[0] == "inference"){
			auto k = 0u; while(k < st.length()-15 && st.substr(k,15) != " burnin-percent") k++;
			if(k < st.length()-15) st = st.substr(0,k);
			if(burnin_frac != BURNIN_FRAC_DEFAULT){
				st += " burnin-percent=" + to_str(burnin_frac);
			}
			lines_raw[j] = st;
			return;
		}
	}
	emsg("Could not find inference");
}


/// Convert number to string
string Output::to_str(double num) const
{
	stringstream ss; ss << num;
	return ss.str();
}


/// Outputs a parameter sample
void Output::param_sample(unsigned int s, unsigned int chain, State &state)
{
	timer[PARAM_OUTPUT] -= clock();
	auto part = state.generate_particle(s,chain,false);

	param_store.push_back(part);
	
	if(sampledir != "" && chain != GEN_PLOT){
		auto nchain = model.details.nchain;
		
		string na = sampledir+"/"+get_file_name("param",chain,nchain,".csv"); 
		
		if(s == 0){
			ofstream fout(na);
			fout << trace_init();		
		}
		
		{
			ofstream fout(na,std::ios::app);
			auto value = param_value_from_vec(part);
			fout << param_output(part,value);
		}
	}
	
	timer[PARAM_OUTPUT] += clock();
}


/// Outputs a parameter sample 
void Output::param_sample(const Particle &part)
{
	timer[PARAM_OUTPUT] -= clock();
	param_store.push_back(part);
	timer[PARAM_OUTPUT] += clock();
}

/// Outputs a parameter sample 
string Output::param_output(const Particle &part, const vector < vector <double> > &value) const
{	
	stringstream ss;
	ss << part.s;

	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		if(par.trace_output){
			if(model.is_symmetric(par)){  // Outputs covariance matrix
				auto L = par.dep[0].list.size();
				for(auto j = 0u; j < L; j++){  // Diagonal elements
					ss << "," << value[th][j*L+j];
				}
	
				for(auto j = 0u; j < L; j++){  // Off-diagonal elements
					for(auto i = j+1; i < L; i++){	
						ss << "," << value[th][j*L+i]; 
					}
				}
			}
			else{			
				for(auto j = 0u; j < par.N; j++){
					if(th >= value.size()) emsg("r");
					if(j >= value[th].size()) emsg("r2");
					ss << "," << value[th][j];
				}
			}
		}
	}
		
	auto dir_out = part.dir_out;
	for(auto i = 0u; i < model.derive.size(); i++){
		const auto &der = model.derive[i];
		if(der.time_dep == false){
			const auto &out = dir_out[i];
			
			if(der.dep.size() == 0){
				ss << "," << out.value_str[0];
			}
			else{
				for(auto j = 0u; j < der.eq.size(); j++){
					ss << "," << out.value_str[j];
				}
			}
		}
	}

	if(model.trans_tree){ 
		const auto &tts = part.trans_tree_stats;
		ss << "," << tts.N_origin << "," << tts.N_inf << "," << tts.N_mut_tree << "," << tts.N_mut_origin << "," << tts.N_unobs << "," << tts.t_root;
	}
	
	for(auto p = 0u; p < model.species.size(); p++){
		if(model.species[p].type == INDIVIDUAL){
			ss << "," << part.species[p].nindividual;
		}
	}
	
	ss << ic_output(part);
	
	const auto &like = part.like;
	
	if(std::isnan(like.init_cond)) emsg("pp");
	ss << "," << like.markov << "," << like.nm_trans << "," << like.ie << "," << like.dist << "," << like.obs << "," << like.genetic_process << "," << like.genetic_obs << "," << like.init_cond;

	if(model.mode == INF) ss << "," << like.prior+like.prior_bounded+like.spline_prior+like.init_cond_prior;
	
	ss << endl;
	return ss.str();
}


/// Outputs a state sample
void Output::state_sample(unsigned int s, unsigned int chain, State &state)
{
	timer[STATE_OUTPUT] -= clock();
	state.check_final_li_wrong();
	auto part = state.generate_particle(s,chain,true);
	state_store.push_back(part);
	timer[STATE_OUTPUT] += clock();
}


/// Outputs a state sample
void Output::state_sample(const Particle &part)
{
	timer[STATE_OUTPUT] -= clock();
	state_store.push_back(part);
	timer[STATE_OUTPUT] += clock();
}

	
/// Outputs a state sample
string Output::state_output(const Particle &part,	vector <string> &ind_key, Hash &hash_ind) const
{
	auto value = param_value_from_vec(part); 
	
	stringstream ss;
	ss << "<<STATE " << part.s << ">>" << endl << endl;
	
	ss << "<PARAMETERS>" << endl;
	
	ss << output_param(value);

	if(use_ind_key == false){
		ss << "<TIMEPOINT " << model.details.t_start << ":" <<  model.details.dt << ":" << model.details.t_end << ">" << endl;
	}
	
	ss << endl;
	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		const auto &ssp = part.species[p];
	
		ss << "<SPECIES \"" << sp.name << "\">" << endl;
		ss << endl;
			
		switch(sp.type){
		case POPULATION:
			{
				ss << "<INITIAL POPULATION>" << endl;
				ss << "compartment,population" << endl;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					ss << sp.comp_gl[c].name << "," << ssp.init_cond_val.cnum[c] << endl;
				}
				ss << endl;
				
				const auto &tp = model.timepoint;
				auto T = tp.size()-1;
				
				if(sp.add_rem_pop_on){
					ss << "<POPULATION CHANGE>" << endl;
					ss << "tp,compartment,change" << endl;
					for(auto ti = 0u; ti < T; ti++){
						for(auto c : sp.add_rem_pop_change[ti]){
							ss << ti << "," << sp.comp_gl[c].name << "," 
									<< sp.add_rem_pop[ti][c] << endl;
						}
					}
					ss << endl;
				}					
				
				ss << "<TRANSITIONS>" << endl;
				// A compressed format
				for(auto i = 0u; i < sp.tra_gl.size(); i++){
					auto pn = replace_arrow(sp.tra_gl[i].name); 
					ss << pn << ":";
					
					auto flag = false;
					auto ti = 0u;
					while(ti < T){
						if(flag == true) ss << ","; 
						flag = true;
						
						if(ssp.trans_num[i][ti] == 0){
							auto num = 0u;
							while(ti < T && ssp.trans_num[i][ti] == 0){ ti++; num++;}
							ss << "Z" << num;
						}
						else{
							ss << ssp.trans_num[i][ti];
							ti++;
						}
					}
					ss << endl;
				}
				ss << endl;
			}
			break;
			
		case INDIVIDUAL:
			{
				ss << "<INDIVIDUALS>" << endl;
				ss << "index,source";
				for(auto i = 0u; i < sp.ind_effect.size(); i++) ss << "," << sp.ind_effect[i].name;
				ss << ",events" << endl;
				for(auto i = 0u; i < ssp.individual.size(); i++){
					const auto &ind = ssp.individual[i];
	
					auto name = ind.name;
					
					if(use_ind_key){
						auto vec = hash_ind.get_vec_string(name);
						auto id = hash_ind.existing(vec);
						if(id == UNSET){
							id = ind_key.size();
							ind_key.push_back(name);
							hash_ind.add(id,vec);
						}
						ss << id << ",";
					}
					else{
						ss << name << ",";
					}
					
					if(i >= sp.nindividual_in && i < sp.nindividual_obs){
						ss << "no";
						for(auto i = 0u; i < sp.ind_effect.size(); i++) ss << "," << ind.exp_ie[i]; 
						ss << ",unobserved";
						
						if(ind.ev.size() != 0) emsg("Should not have an event");
					}
					else{
						if(ind.ev.size() == 0) emsg("Should have an event");
						
						const auto &ev = ind.ev[0];
						switch(ev.type){
						case ENTER_EV: ss << "no"; break;
						case NM_TRANS_EV: case M_TRANS_EV:
							{
								const auto &tr = sp.tra_gl[ev.tr_gl];
								if(tr.variety != SOURCE_TRANS) emsg("Should be source");
								ss << "yes";
							}
							break;
							
						case LEAVE_EV: case MOVE_EV: 
							emsg("Should not start with this event");
							break;
						}
				
						for(auto i = 0u; i < sp.ind_effect.size(); i++) ss << "," << ind.exp_ie[i]; 
						ss << ",";
						auto c = UNSET;
						for(auto e = 0u; e < ind.ev.size(); e++){
							const auto &ev = ind.ev[e];
							auto cnew = c;
							
							if(e != 0) ss << " ";
							switch(ev.type){
							case ENTER_EV:
								{
									if(c != UNSET) emsg("Inconsistent");
									cnew = ev.c_after; 
									ss << sp.comp_gl[cnew].name;
									infection_info(ev.ind_inf_from,p,part,ss);
								}
								break;

							case LEAVE_EV: 
								cnew = UNSET; 
								ss << "-"; 
								break;
								
							case MOVE_EV: 
								if(c != UNSET){
									auto cl = ev.cl;
									cnew = sp.update_c_comp(c,cl,ev.move_c);
									ss << "move->" << sp.cla[cl].comp[ev.move_c].name;
									infection_info(ev.ind_inf_from,p,part,ss);
								}
								break;
								
							case M_TRANS_EV: case NM_TRANS_EV: 
								{
									const auto &trg = sp.tra_gl[ev.tr_gl];
									if(trg.i == UNSET) ss << replace_arrow(trg.name);
									else{
										const auto &tr = sp.cla[trg.cl].tra[trg.tr];
										ss << replace_arrow(tr.name);
									}
									
									infection_info(ev.ind_inf_from,p,part,ss);
									
									if(c != trg.i) emsg("Inconsist");
									cnew = trg.f;
								}
								break;
							}
							ss << ":" << model.calc_t(ev.tdiv);
							c = cnew;
						}
					}
					ss << endl;
				}
				ss << endl;
			}
			break;
			
		default:
			emsg("op error");
			break;
		}
		
		if(cum_diag && model.mode == INF){
			ss << "<TRANSDISTPROB>" << endl;
				
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				ss << replace_arrow(sp.tra_gl[tr].name) << ":";
				for(auto i = 0u; i < H_BIN; i++){
					if(i != 0) ss << "|";
					ss << ssp.cum_prob_dist[tr][i];
				}
				ss << endl;
			}
			ss << endl;
		}
	}
	
	auto j = 0u; while(j < model.derive.size() && model.derive[j].time_dep == false) j++;
	if(j < model.derive.size()){
		auto dir_out = part.dir_out;
		
		ss << "<DERIVED>" << endl;
		for(auto i = 0u; i < model.derive.size(); i++){
			const auto &der = model.derive[i];
			const auto &out = dir_out[i];
			
			if(der.time_dep == true){
				ss << der.name;
				if(der.dep.size() == 0){
					ss << "," << out.value_str[0] << endl;
				}
				else{
					ss << endl;
					ss << "  ";
					for(auto k = 0u; k < der.dep.size(); k++){
						if(k != 0) ss << ",";
						ss << der.dep[k].index_with_prime;
					}
					ss << ",Value" << endl;
						
					for(auto j = 0u; j < der.eq.size(); j++){
						ss << "  ";
						for(auto k = 0u; k < der.dep.size(); k++){
							const auto &dp = der.dep[k];
							auto m = (unsigned int)(j/dp.mult)%dp.list.size();
							if(k != 0) ss << ",";
							ss << dp.list[m];
						}
						ss << "," << out.value_str[j] << endl;
					}
				}
			}
		}
	}
	
	if(model.trans_tree){                   // Outputs the phylogenetic tree
		ss << "<TRANSTREE>" << endl;
		//if(part.inf_node.size() == 0) emsg(" zero");
		for(const auto &in : part.inf_node){
			ss << part.species[in.p].individual[in.i].name  << "|" << model.calc_t(in.tdiv_start) << "|";	
			
			auto n_from = in.from.node;
			if(n_from == ENTER_INF || n_from == OUTSIDE_INF){
				if(n_from == ENTER_INF) ss << "ENT|";
				else ss << "OUT|";
				auto num = part.inf_origin[in.from.index].mut_num;
				if(num == UNSET) ss << "NA";
				else ss << num;
			}
			else{
				ss << n_from << "|" << in.from.index;
			}
			
			for(const auto &iev : in.inf_ev){
				ss << ",";
				switch(iev.type){
				case INFECT_OTHER: ss << "INF"; break;
				case GENETIC_OBS: ss << "OBS"; break;
				}
				
				ss << "|" << iev.index << "|" << model.calc_t(iev.tdiv) << "|";
				if(iev.mut_num == UNSET) ss << "NA";
				else ss << iev.mut_num;
			}
			ss << endl;
		}
		ss << endl;
	}
	
	return ss.str();
}


/// Generate test to indicate an infection
void Output::infection_info(const IndInfFrom &iif, unsigned int p, const Particle &part, stringstream &ss) const
{
	if(iif.p != UNSET){
		ss << "[";
		switch(iif.p){
		case ENTER_INF: ss << "ENT_INF"; break; 
		case OUTSIDE_INF: ss << "OUT_INF"; break;
		default:
			{
				auto pp = iif.p;
				if(pp != p) ss << model.species[pp].name << "|";
				ss << part.species[pp].individual[iif.i].name;
			}
			break;
		}
		ss << "]";
	}
}


/// Generates the final output files
string Output::generate_state_head(const vector <string> &ind_key) const
{
	stringstream ss;
	
	if(use_ind_key){	
		ss << "{" << endl;
		ss << "  # Time divisions" << endl;
		ss << "  timepoint " << model.details.t_start << ":" <<  model.details.dt << ":" << model.details.t_end << endl;
		
		auto out_fl = false;
		for(auto &par : model.param){ if(par.spline_outside) out_fl = true;}
		if(out_fl == true){
			ss << "  # Spline out of range" << endl;
			for(auto &par : model.param){
				if(par.spline_outside){
					ss << "  spline-out " << par.name << " ";
					const auto &dep = par.dep;
					const auto &list = dep[dep.size()-1].list;
					for(auto k = 0u; k < list.size(); k++){
						if(k != 0) ss << ",";
						ss << list[k];
					}
					ss << endl;
				}
			}
		}
		
		if(ind_key.size() > 0){
			ss << "  # Individual index" << endl;
			for(auto i = 0u; i < ind_key.size(); i++){
				ss << "  " << i << ":" << ind_key[i] << endl;
			}
		}
		ss << "}" << endl << endl;
	}
	return ss.str();
}


/// Generates values from parameter vector
vector < vector <double> > Output::param_value_from_vec(const Particle &pa) const 
{
	auto param_val = model.get_param_val(pa);
	//cout << "PRINT" << endl;
	//model.print_param(param_val);
	
	model.add_tvreparam(param_val,pa.param_val_tvreparam);
	
	const auto &val = param_val.value;
	vector < vector <double> > value;
	
	value.resize(model.param.size());
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		if(par.trace_output){
			value[th].resize(par.N,UNSET);
			if(par.variety != CONST_PARAM){
				for(auto j = 0u; j < par.N; j++){
					const auto &er = par.element_ref[j];
					auto ind = er.index;
					if(ind != UNSET){
						if(er.cons){
							value[th][j] = par.constant.value[ind];
						}
						else{
							const auto &ele = par.element[ind];
							auto k = ele.param_vec_ref;
							if(k != UNSET) value[th][j] = val[k];
						}
					}
				}
			}
		}
	}
	
	if(val.size() != model.param_vec.size()) emsg("param_val the wrong size");
	
	/*
	for(auto i = 0u; i < model.param_vec.size(); i++){
		const auto &mpv = model.param_vec[i];
		auto th = mpv.th;
		if(model.param[th].trace_output){
			value[mpv.th][mpv.index] = val[i]; 
		}
	}
	*/
	
	if(false){
		for(auto th = 0u; th < model.param.size(); th++){
			const auto &par = model.param[th];
			cout << par.name << endl;
			for(auto j = 0u; j < par.N; j++){
				cout << value[th][j]  << ",";
			}
			cout << " va" << endl;
		}
	}
	
	const auto &precalc = param_val.precalc;
	
	// Fills in any unset values
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		if(par.trace_output){
			for(auto j = 0u; j < par.N; j++){
				if(value[th][j] == UNSET){	
					auto val = par.get_value(j);
					if(val != UNSET){
						value[th][j] = val;
					}
					else{
						auto eq = par.get_eq_ref(j);
						if(eq != UNSET){			
							value[th][j] = model.eqn[eq].calculate_param(precalc);
						}
					}
				}
			}
		}
	}
	
	return value;
}


/// Generates parameters values from those already stored in the parameter
vector < vector <double> > Output::param_value_from_value() const 
{
	vector < vector <double> > value;
	
	value.resize(model.param.size());
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		if(par.variety == CONST_PARAM){
			value[th].resize(par.N);
			for(auto i = 0u; i < par.N; i++) value[th][i] = par.get_value(i);
		}
	}

	return value;
}


/// Outputs a string version of parameters which are either constant or non-constant 
string Output::output_param(const vector < vector <double> > &value) const
{
	stringstream ss;
	
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		
		if(par.trace_output){
			auto pn = replace_arrow(par.name);
			pn = add_escape_char(pn);
			ss << "\"" << pn << "\"";
			if(par.dep.size() == 0){
				ss << "," << value[th][0] << endl;
			}
			else{
				ss << endl;
				ss << "  ";
				if(model.is_symmetric(par)){       // Outputs covariance matrix
					auto L = par.dep[0].list.size();
					for(auto i = 0u; i < L; i++){
						if(i != 0) ss << ",";
						ss << par.dep[0].list[i];
					}
					ss << endl;
					
					for(auto j = 0u; j < L; j++){ 
						ss << "  ";
						for(auto i = 0u; i < L; i++){
							if(i != 0) ss << ",";
							if(i < j) ss << ".";
							else ss << value[th][j*L+i]; 
						}
						ss << endl;
					}
				}
				else{
					for(auto k = 0u; k < par.dep.size(); k++){
						if(k != 0) ss << ",";
						ss << par.dep[k].index_with_prime;
					}
					ss << ",Value" << endl;
						
					for(auto j = 0u; j < value[th].size(); j++){
						ss << "  ";
						for(auto k = 0u; k < par.dep.size(); k++){
							const auto &dp = par.dep[k];
							auto m = (unsigned int)(j/dp.mult)%dp.list.size();
							if(k != 0) ss << ",";
							ss << dp.list[m];
						}
						ss << "," << value[th][j] << endl;
					}
				}
			}
		}
	}
	/*
	const auto &like = part.like;
	ss << "L^markov," << like.markov << endl;
	ss << "L^non-markov," << like.nm_trans << endl;
	ss << "L^ie," << like.ie << endl;
	ss << "L^dist," << like.dist << endl;
	ss << "L^obs," << like.obs << endl;
	ss << "L^init," << like.init_cond << endl;
	ss << "Prior," << like.prior << endl;
	*/
	
	ss << endl;
	
	return ss.str();
}


/// Gets the particles for a given chain and deletes from list
vector <Particle> Output::get_part_chain(unsigned int chain, vector <Particle> &part) const
{
	vector <Particle> part_ch;
	
	vector <unsigned int> hole;
	auto j = 0u;
	
	for(auto i = 0u; i < part.size(); i++){
		auto &pa = part[i];
		if(pa.chain == chain){
			part_ch.push_back(pa);
			clear_part(pa);
			hole.push_back(i);
		}
		else{
			if(j < hole.size()){
				part[hole[j]] = pa;
				clear_part(pa);
				j++;
			}
		}
	}
	
	if(j < hole.size()){
		part.resize(part.size()-(hole.size()-j));
	}
	
	return part_ch;
}


/// Clears a particle
void Output::clear_part(Particle &pa) const
{
	pa.param_val_prop.clear();
	pa.species.clear();
	pa.dir_out.clear();
	pa.inf_origin.clear();
	pa.inf_node.clear();
}


/// Prints out observation events 
string Output::print_obs_data(unsigned int p, const vector <ObsData> &obs) const 
{
	stringstream ss;
	
	const auto &sp = model.species[p];
	
	ss << "Observations:" << endl;
	for(auto &ob : obs){
		const auto &claa =  sp.cla[ob.cl];
		
		ss << model.calc_t(ob.tdiv) << " ";
		switch(ob.type){
		case OBS_TRANS_EV: 
			ss << "OBS_TRANS_EV " << sp.obs_trans[ob.ref].name;
			break;

		case OBS_SOURCE_EV: 
			ss << "OBS_SOURCE_EV " << sp.obs_trans[ob.ref].name;
			break;

		case OBS_SINK_EV: 
			ss << "OBS_SINK_EV " << sp.obs_trans[ob.ref].name;
			break;
			
		case OBS_COMP_EV:
			{
				ss << "OBS_COMP_EV ";
				const auto &comp = claa.comp;
				if(ob.c_exact != UNSET){
					ss << comp[ob.c_exact].name;
				}
				else{				
					for(auto c = 0u; c < comp.size(); c++){
						auto eq = sp.obs_eqn[ob.obs_eqn_ref[c]];
						ss << comp[c].name << ":" << model.eqn[eq].te_raw;
					}
				}
			}
			break;
			
		case OBS_TEST_EV:
			{
				ss << "OBS_TEST_EV res=";
				if(ob.test_res == true) ss << "+ve"; else ss << "-ve";
				//const auto &ds = sp.source[ob.so];
				
				
				ss << "  Se=" << ob.Se_eqn.te << "  Sp=" << ob.Sp_eqn.te;
			}
			break;
			
		default: emsg("event data not right:"+ob.type); break;
		}
		ss << "   ";
	}
	ss << endl << endl;
	
	return ss.str();
}


/// Prints information about individuals
void Output::print_individuals(unsigned int N, unsigned int p, const State &state) const
{
	const auto &sp = model.species[p];
	const auto &ssp = state.species[p];
	
	if(N == UNSET) N = sp.individual.size();
	
	for(auto i = 0u; i < N; i++){
		const auto &ind = sp.individual[i];
		
 		cout << endl << ind.name << " Individual" << endl;
		ssp.print_event("Ev",ssp.individual[i]);
		state.print_ev_data("Event Data",ind.ev,p);
		cout << print_obs_data(p,ind.obs);
	}
}


/// Sets diagnostics
void Output::set_diagnostics(unsigned int ch, string diag)
{
	Diagnostic di; di.ch = ch; di.te = diag;
	diagnostic_store.push_back(di);
}


/// Generates the final outputs
void Output::end(string file, unsigned int total_cpu)
{
	percentage_start(OUTPUT_PER);
		
	if(com_op) cout << "<OUTPUTTING>" << endl;
	
	if(com_op == true){
		cout << "<<OUTPUT FILE>>" << endl;
		if(mpi.core == 0){
			cout << "# PROCESSED USING BICI " << bici_version << endl;
			for(auto li : lines_raw){
				cout << li << endl;
			}
		}
	}
	
	ofstream fout;
	
	if(com_op == false && op()){
		fout.open(file);

		fout << "# PROCESSED USING BICI " << bici_version << endl;
		for(auto li : lines_raw){
			fout << li << endl;
		}
	}
	
	auto nchain = model.details.nchain;

	vector < vector < vector < vector <double> > > > param_samp;
	param_samp.resize(nchain);

	for(auto ch = 0u; ch < nchain; ch++){
		auto frac = double(ch)/nchain;
		percentage(frac*25,100);
		
		auto part = get_part_chain(ch,param_store);
		
#ifdef USE_MPI	
		mpi.transfer_particle(part);
#endif
	
		if(op() && part.size() > 0){
			number_part(part);
			output_trace(ch,part,param_samp,fout);
		}
	}

	vector <string> final_warning;
	
	if(op() && model.mode == INF && com_op == false){
		output_param_statistics(param_samp,fout,final_warning);
	}
	
	auto trans_diag = trans_diag_init(); 
	
	auto trans_diag_on = false; if(model.mode == INF) trans_diag_on = true;
	
	for(auto ch = 0u; ch < model.details.nchain; ch++){
		auto frac = double(ch)/nchain;
		percentage(25+frac*25,100);

		auto part = get_part_chain(ch,state_store);
		
#ifdef USE_MPI	
		mpi.transfer_particle(part);
#endif
	
		if(op() && part.size() > 0){
			number_part(part);
			if(trans_diag_on) trans_diag_add(trans_diag,part); 

			output_state(ch,part,fout);
		}
	}

	if(op() && trans_diag_on){
		output_trans_diag(trans_diag,fout);
	}
	
	auto alg = model.details.algorithm;
	if(alg == PAS_MCMC || alg == ABC_SMC_ALG){
		auto part = get_part_chain(GEN_PLOT,param_store);
		
#ifdef USE_MPI	
		mpi.transfer_particle(part);
#endif
		
		if(op() && part.size() > 0){
			output_generation(part,fout);
		}
	}

	auto alg_warn_flag = false;
	
	if(model.details.diagnostics_on){
		auto diagnostic = diagnostic_store;
	
#ifdef USE_MPI	
		mpi.transfer_diagnostic(diagnostic);
#endif

		if(op() && diagnostic.size() > 0){
			output_diagnostic(diagnostic,alg_warn_flag,fout);
		}
	}

	output_add_ind_warning(final_warning);

	output_rate_warning(total_cpu,50,100,final_warning);
	
	output_spline_out_warning(final_warning);

	for(const auto &der : model.derive){
		auto st = der.func.warn; 
		if(der.func.on && st != "") final_warning.push_back(st);
	}

	if(mpi.core == 0){ // Warnings for data out of range
		for(const auto &sp : model.species){
			for(auto warn : sp.data_warning){
				auto wa = warn;
				if(model.species.size() > 1) wa += "For species '"+sp.name+"': "+wa;
				final_warning.push_back(wa);
			}
		}
	}
	
	if(op()){
		for(auto te : final_warning) add_warning(te,fout);
		if(alg_warn_flag) add_warning("This run has generated algorithm warnings. Please check the diagnostic file(s) for details",fout);
	}
	
	if(com_op == true) cout << "<<END>>" << endl;
}


/// Outputs trace 
void Output::output_trace(unsigned int ch, const vector <Particle> &part, vector < vector < vector < vector <double> > > > &param_samp, ofstream &fout) const
{
	auto nchain = model.details.nchain;
	
	auto burn = 0u;
	const auto &de = model.details;
	switch(de.algorithm){
	case DA_MCMC: case PAS_MCMC: burn = double(de.sample*de.burnin_frac)/100; break;
	default: break;
	}

	string param_out;
	
	param_out += trace_init();
	for(const auto &pa : part){
		auto value = param_value_from_vec(pa);

		if(model.mode == INF && pa.s >= burn){
			param_samp[ch].push_back(value);				
		}
		param_out += param_output(pa,value);
	}
	
	string param_out_file;

	if(sampledir != ""){
		param_out_file = get_file_name("param",ch,nchain,".csv"); 

		ofstream pout(sampledir+"/"+param_out_file);
		check_open(pout,param_out_file);
		pout << param_out;

		param_out_file = sampledir_rel+"/"+param_out_file;
	}
	else{  // Embeds output into file
		param_out_file = "[["+endli+param_out+"]]";
	}

	stringstream ss;
	switch(model.mode){
	case SIM:
		ss << "param-sim file=\"" << param_out_file << "\"" << endl;
		break;

	case INF:
		ss << "param-inf chain=" << ch << " file=\""+param_out_file+"\"" << endl;
		break;
	
	case PPC:
		ss << "param-post-sim file=\"" << param_out_file << "\"" <<  endl;
		break;
	
	default: emsg("Should not be default12"); break;
	}
	ss << endl;
	
	if(com_op == true) cout << ss.str();
	else fout << ss.str();
}


/// Gets a file name based on a root and chain number 
string Output::get_file_name(string root, unsigned int ch, unsigned int nchain, string end) const
{
	auto name = root;
	if(ch != UNSET && ch != GEN_PLOT && nchain > 1 && model.mode != PPC) name += "_"+tstr(ch);
	name += end;
	
	return name;
}


/// Calculate the effective sample size for a list of numbers
string Output::get_effective_sample_size(vector <double> vec) const
{
	auto N = vec.size();
	
	double min = LARGE, max = -LARGE;
	for(auto i = 0u; i < N; i++) {
		auto val = vec[i];
		if(val < min) min = val;
		if(val > max) max = val;
	}
	if(min == max || N < 2) return "NA";
		
	auto av = 0.0, av2 = 0.0;
	for(auto s = 0u; s < N; s++){
		auto val = vec[s]; av += val; av2 += val*val;
	}
	auto num = av2/N - (av/N)*(av/N); if(num < 0) num = 0;
	auto mean = av/N, sd = sqrt(num);
	
	for(auto s = 0u; s < N; s++) vec[s] = (vec[s]-mean)/sd;
	
	auto sum = 1.0;
	for(auto d = 1u; d < N/2; d++){
		auto a = 0.0; for(auto s = 0u; s < N-d; s++) a += vec[s]*vec[s+d]; 
		auto cor = a/(N-d); if(cor < 0) break;
		sum += 2*cor;			
	}
	
	return tstr((unsigned int)(N/sum));
}	


/// Returns the Gelman-Rubin statistic
string Output::get_Gelman_Rubin_statistic(const vector < vector <double> > &cha) const
{
	auto C = cha.size();
	if(C == 0) return "no chain";
	if(C == 1) return "NA";
	
	auto N = cha[0].size();
	if(N == 0) return "no data";
		
	vector <double> mu(C), vari(C);
	
	{ // Checks to see if values along chain are all equal
		for(auto ch = 0u; ch < C; ch++){ 
			auto val_basic = cha[ch][0];
			auto i = 0u; while(i < N && cha[ch][i] == val_basic) i++;
			if(i == N) return "NA";
		}
	}
	
	auto muav = 0.0;
	for(auto ch = 0u; ch < C; ch++){ 
		auto valav = 0.0; for(auto i = 0u; i < N; i++) valav += cha[ch][i]/N;
		auto varr = 0.0; for(auto i = 0u; i < N; i++) varr += (cha[ch][i]-valav)*(cha[ch][i]-valav)/(N-1);
		mu[ch] = valav;
		vari[ch] = varr;
		muav += mu[ch]/C;
	}
	
	auto W = 0.0; for(auto ch = 0u; ch < C; ch++) W += vari[ch]/C;
	
	auto B = 0.0; for(auto ch = 0u; ch < C; ch++) B += (mu[ch]-muav)*(mu[ch]-muav)*N/(C-1);
	
	return tstr(sqrt(((1-1.0/N)*W + B/N)/W));
}


/// Generates parameter statistics
string Output::param_stat(unsigned int th, unsigned int i, const vector < vector < vector < vector <double> > > > &param_samp, vector <Warn> &ESS_warn, vector <Warn> &GR_warn) const
{
	auto nchain = param_samp.size();
	
	vector <double> vec;
	vector < vector <double> > cha;
	cha.resize(nchain);
	
	for(auto ch = 0u; ch < nchain; ch++){
		const auto &psch = param_samp[ch];
		for(auto s = 0u; s < psch.size(); s++){
			auto val = psch[s][th][i];
			vec.push_back(val);
			cha[ch].push_back(val);
		}
	}
	
	auto ESS = get_effective_sample_size(vec);
	auto GR = get_Gelman_Rubin_statistic(cha);
	
	auto ESS_num = number(ESS);
	if(ESS_num != UNSET && ESS_num < ESS_THRESH){
		if(find_in(ESS_warn,th) == UNSET){
			Warn wa; wa.th = th; wa.num = ESS_num;
			ESS_warn.push_back(wa);
		}
	}			
	
	auto GR_num = number(GR);
	if(GR_num != UNSET && GR_num > GR_THRESH){
		if(find_in(GR_warn,th) == UNSET){
			Warn wa; wa.th = th; wa.num = GR_num;
			GR_warn.push_back(wa);
		}
	}			
	
	auto stat = get_statistic(vec);
	
	stringstream ss;
	ss << "," << stat.mean << "," << stat.sd << "," << stat.CImin << "," << stat.CImax << "," << ESS << "," << GR;

	return ss.str();
}


/// Outputs a table of parameter statistics
void Output::output_param_statistics(const vector < vector < vector < vector <double> > > > &param_samp, ofstream &fout, vector <string> &final_warning) const
{
	stringstream sss;
	
	sss << "name,mean,sd,CI min,CI max,ESS,GR" << endl;
	
	vector <Warn> ESS_warn, GR_warn;
	
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		
		if(par.trace_output){
			auto pn = add_escape_char(par.name);
			
			if(par.dep.size() == 0){	
				sss << "\"" << pn << "\"" << param_stat(th,0,param_samp,ESS_warn,GR_warn) << endl;
			}
			else{
				for(auto j = 0u; j < par.N; j++){
					sss << "\"" << pn << "_";
					for(auto k = 0u; k < par.dep.size(); k++){
						const auto &dp = par.dep[k];
						auto m = (unsigned int)(j/dp.mult)%dp.list.size();
						if(k != 0) sss << ",";
						sss << dp.list[m];
					}
					sss << "\"";
					sss << param_stat(th,j,param_samp,ESS_warn,GR_warn) << endl;
				}
			}
		}
	}
	
	auto stat_table = sss.str();

	string stat_out_file;
	
	if(sampledir != ""){
		stat_out_file = "param-stats-inf.csv";

		ofstream pout(sampledir+"/"+stat_out_file);
		check_open(pout,stat_out_file);
		pout << stat_table;

		stat_out_file = sampledir_rel+"/"+stat_out_file;
	}
	else{  // Embeds output into file
		stat_out_file = "[["+endli+stat_table+"]]";
	}
	
	stringstream ss;
	ss << "param-stats-inf file=\""+stat_out_file+"\"" << endl;
	ss << endl;
	
	if(com_op == true) cout << ss.str();
	else fout << ss.str();
	
	if(ESS_warn.size() > 0 || GR_warn.size() > 0){
		stringstream ss;
		
		ss << "MCMC diagnostics suggests it has not been run for long enough. ";
		if(ESS_warn.size() > 0){
			auto ESS_min = LARGE;
			for(const auto &wa : ESS_warn){ if(wa.num < ESS_min) ESS_min = wa.num;}
			
			auto sa = (unsigned int)(1.1*model.details.sample*double(ESS_THRESH)/ESS_min);
			ss << "An estimated " << sa << " updates are required. ";
		}
		
		//ss << endl;
		
		if(ESS_warn.size() > 0){
			string pte;
			for(const auto &wa : ESS_warn){
				if(pte != "") pte += ", ";
				pte += model.param[wa.th].full_name;
			}
			
			//ss << endl;
			if(ESS_warn.size() == 1) ss << "Parameter";
			else ss << "Parameters";
			ss << " below the " << ESS_THRESH << " ESS threshold: " << pte << ". ";
		}
		
		if(GR_warn.size() > 0){
			string pte;
			for(const auto &wa : GR_warn){
				if(pte != "") pte += ", ";
				pte += model.param[wa.th].full_name;
			}
			
			//ss << endl;
			if(GR_warn.size() == 1) ss << "Parameter";
			else ss << "Parameters";
		
			ss << " above the " << GR_THRESH << " GR threshold: " << pte << ". ";
		}
		
		final_warning.push_back(trim(ss.str()));
	}
}
	

/// Outputs state samples
void Output::output_state(unsigned int ch, const vector <Particle> &part, ofstream &fout) const
{
	auto nchain = model.details.nchain;
	
	string state_out;

	vector <string> ind_key;
	Hash hash_ind;
	for(const auto &pa : part) state_out += state_output(pa,ind_key,hash_ind);
	
	auto state_head = generate_state_head(ind_key);
	
	string state_out_file;

	if(sampledir != ""){
		state_out_file = get_file_name("state",ch,nchain,".txt"); 
		
		ofstream pout(sampledir+"/"+state_out_file);
		check_open(pout,state_out_file);
		pout << state_head << state_out;

		state_out_file = sampledir_rel+"/"+state_out_file;
	}
	else{  // Embeds output into file
		state_out_file = "[["+endli+state_head+state_out+"]]";
	}

	stringstream ss;
	switch(model.mode){
	case SIM:
		ss << "state-sim file=\"" << state_out_file << "\"" << endl;
		break;

	case INF:
		ss << "state-inf chain=" << ch << " file=\""+state_out_file+"\"" << endl;
		break;
	
	case PPC:
		ss << "state-post-sim file=\"" << state_out_file << "\"" <<  endl;
		break;
	
	default: emsg("Should not be default12"); break;
	}
	ss << endl;
	
	if(com_op == true) cout << ss.str();
	else fout << ss.str();
}

	
/// Outputs results for trans_diag
void Output::output_trans_diag(const vector <TransDiagSpecies> &trans_diag, ofstream &fout) const
{
	stringstream sstd;
	
	auto T = model.details.T;
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		if(sp.tra_gl.size() > 0){
			if(model.nspecies > 1) sstd << "<SPECIES \"" << sp.name << "\">" << endl;
			const auto &exp_num = trans_diag[p].exp_num;
			auto n = trans_diag[p].n;
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				sstd << replace_arrow(sp.tra_gl[tr].name)  << ":";
				for(auto ti = 0u; ti < T; ti++){
					if(ti != 0) sstd << ",";
					sstd << exp_num[tr][ti]/n;
				}
				sstd << endl;
			}
		}
	}
	
	auto trans_diag_out = sstd.str();
	string trans_diag_out_file;
	
	if(sampledir != ""){
		trans_diag_out_file = "diagnostic/trans_diag.txt";

		ofstream pout(sampledir+"/"+trans_diag_out_file);
		check_open(pout,trans_diag_out_file);
		pout << trans_diag_out;

		trans_diag_out_file = sampledir_rel+"/"+trans_diag_out_file;
	}
	else{  // Embeds output into file
		trans_diag_out_file = "[["+endli+trans_diag_out+"]]";
	}

	stringstream ss;
	ss << "trans-diag-inf file=\"" << trans_diag_out_file << "\"" << endl;
	ss << endl;
	
	if(com_op == true) cout << ss.str();
	else fout << ss.str();
}


/// Adds an output warning
void Output::add_warning(string err_msg, ofstream &fout) const
{
	string line;
	
	switch(model.mode){
	case SIM: line = "warning-sim"; break;
	case INF: line = "warning-inf"; break;
	case PPC: line = "warning-post-sim"; break;
	case MODE_UNSET: emsg("op problem"); break;
	}
	
	line += " text=\"[["+endli;
	line += err_msg+endli;
	line += "]]\""+endli+endli;
	
	if(com_op == true){
		cout << line;
	}
	else{
		cout << endl;
		display_warning(err_msg);
		cout << endl;
		fout << line;
	}
}

	
/// Outputs a warning if individuals are in the system but not explicitly added
void Output::output_add_ind_warning(vector <string> &final_warning) const
{
	for(const auto &sp : model.species){
		if(sp.type == INDIVIDUAL){
			string te = "";
			auto added = false;
			for(auto &ind : sp.individual){
				if(ind.ev.size() != 0){
					if(ind.ev[0].type == ENTER_EV) added = true;
					else{
						if(te.length() > 50){ te += "..."; break;}	
						if(te != "") te += ", ";
						te += "'"+ind.name+"'";
					}
				}					
			}
			
			if(added == true && te != ""){
				auto msg = "The following individuals are not explicitly added to the system (is this correct?): "+te;
				final_warning.push_back(msg);
			}
		}
	}
}
		
			
/// Works out if time-step is too small
void Output::output_rate_warning(unsigned int total_cpu, unsigned int per_start, unsigned int per_end, vector <string> &final_warning) const 
{
	State state(model);
	state.init();
	
	vector < vector < vector <double> > > rate_sum;
	double nrate_sum = 0;
	
	auto T = state.T;
	auto dt = model.details.dt;
	
	rate_sum.resize(model.nspecies); 
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];	
		rate_sum[p].resize(sp.tra_gl.size()); 
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			rate_sum[p][tr].resize(T,0); 
		}
	}
	
	// Randomly selects around 10 state samples on which estimate are made
	auto samp_tot = state_store.size()*mpi.ncore;
	auto step = samp_tot/10;
	if(step == 0) step = 1;
	
	auto rmax = 0.0;
	for(auto k = 0u; k < state_store.size(); k += step){
		auto frac = double(k)/state_store.size();
		percentage(per_start+frac*(per_end-per_start),100);
		
		const auto &pa = state_store[k];
	
		state.set_particle(pa,false);
		for(auto p = 0u; p < model.nspecies; p++){
			const auto &sp = model.species[p];
			auto rate = state.get_population_rates(p);
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				for(auto ti = 0u; ti < T; ti++){
					rate_sum[p][tr][ti] += rate[tr][ti];
				}
			}
		}
		nrate_sum++;
	}
	

#ifdef USE_MPI	
	for(auto p = 0u; p < model.nspecies; p++){
		mpi.sum(rate_sum[p]);
	}
	mpi.sum(nrate_sum);
#endif

	percentage_end();

	if(op()){
		auto rate_thresh = RATE_RECOMMEND/dt;
		
		string err_msg = "";
		for(auto p = 0u; p < model.nspecies; p++){		
			const auto &sp = model.species[p];
	
			auto K = sp.ncla;
		
			vector < vector < vector <unsigned int> > > tra_good, tra_bad;
			tra_good.resize(K); tra_bad.resize(K); 
			for(auto cl = 0u; cl < K; cl++){
				const auto &claa = sp.cla[cl];
				tra_good[cl].resize(claa.ntra); tra_bad[cl].resize(claa.ntra); 
			}

			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				const auto &tra = sp.tra_gl[tr];
				
				auto fl = false;
				if(tra.i == UNSET){ // Source transition
					
				}
				else{
					for(auto ti = 0u; ti < T; ti++){
						auto val = rate_sum[p][tr][ti]/nrate_sum;
						if(val > rmax) rmax = val;
						if(val > rate_thresh) fl = true;
					}
				}
				if(fl == false) tra_good[tra.cl][tra.tr].push_back(tr);
				else tra_bad[tra.cl][tra.tr].push_back(tr);
			}
			
			string err = "";
			for(auto cl = 0u; cl < K; cl++){
				const auto &claa = sp.cla[cl];
				for(auto tr = 0u; tr < claa.ntra; tr++){
					const auto &tbad = tra_bad[cl][tr];
					const auto &tgood = tra_good[cl][tr];
					if(tbad.size() > 0){
						if(err != "") err += ", ";
						err += claa.tra[tr].name;
						if(tgood.size() != 0){
							string filt = "";
							for(auto k = 0u; k < tbad.size(); k++){
								const auto &trg = sp.tra_gl[tbad[k]];
								auto c = trg.i; if(c == UNSET) c = trg.f;
								const auto &co = sp.comp_gl[c];
								
								string pst = "";
								for(auto cl2 = 0u; cl2 < K; cl2++){
									if(cl2 != cl){
										if(pst != "") pst += "|";
										pst += sp.cla[cl2].comp[co.cla_comp[cl2]].name;
									}
								}
								if(filt != "") filt += ", ";
								filt += pst;
							}
							filt = trunc(filt,40);
							err += " (for "+filt+")";
						}
					}
				}
			}
			
			if(err != ""){
				if(model.nspecies > 1) err_msg += "For species '"+sp.name+"': ";
				err_msg += err+". ";
			}
		}
	
		if(err_msg != ""){
			err_msg = "High transition rates mean that there is a potential for a finite time-step discretisation error. It is recommended that this model is run with a time-step below '"+tstr(RATE_RECOMMEND/rmax,2)+"'. The following transitions are affected: "+err_msg;
		}
		else{
			if(total_cpu >= 100 && rmax != 0 && rmax < 0.1*rate_thresh){
				auto tmax = model.details.t_end-model.details.t_start;
				auto dt = RATE_RECOMMEND/rmax;
				if(dt > tmax/10) dt = tmax/10;
				
				err_msg = "This is being run with a relatively small time-step. Analysis suggests it could be reliably run up to a time-step '"+tstr(dt,2)+"'.";
			}
		}
		
		if(err_msg != ""){
			err_msg = trunc(trim(err_msg),300);
			
			if(final_warning.size()== 0){		
				final_warning.push_back(err_msg);
			}
		}
	}
}


/// Outputs a warning if a spline is outside time range
void Output::output_spline_out_warning(vector <string> &final_warning) const
{
	string te = "";
	for(const auto &par : model.param){
		if(par.spline_outside_par){
			if(te != "") te += ", ";
			te += par.full_name;
		}
	}
	
	if(te != ""){
		final_warning.push_back("The following splines were truncated because they were outside the time range: "+te);
	}
}


/// Outputs information about generations (for PAS-MBP and ABC-SMC)
void Output::output_generation(const vector <Particle> &part, ofstream &fout) const
{
	string param_out;
			
	string content;
	for(const auto &pa : part){
		auto value = param_value_from_vec(pa);
		content += param_output(pa,value);
	}
	param_out += generation_average(trace_init(),content);

	string param_out_file;

	if(sampledir != ""){
		param_out_file = "diagnostic/generation.csv";

		ofstream pout(sampledir+"/"+param_out_file);
		check_open(pout,param_out_file);
		pout << param_out;

		param_out_file = sampledir_rel+"/"+param_out_file;
	}
	else{  // Embeds output into file
		param_out_file = "[["+endli+param_out+"]]";
	}

	stringstream ss;
	ss << "generation-inf file=\""+param_out_file+"\"" << endl;
	ss << endl;
	
	if(com_op == true) cout << ss.str();
	else fout << ss.str();
}


/// Outputs diagnostic information
void Output::output_diagnostic(const vector <Diagnostic> &diagnostic, bool &alg_warn_flag, ofstream &fout) const
{
	for(const auto &di :  diagnostic){
		if(begin_str(di.te,ALG_WARN)) alg_warn_flag = true;
		
		string diag_out_file;
	
		if(diagdir != ""){
			auto nchain = model.details.nchain;
			
			diag_out_file = get_file_name("info",di.ch,nchain,".txt");  

			ofstream pout(diagdir+"/"+diag_out_file);
			check_open(pout,diag_out_file);
			pout << di.te;

			diag_out_file = sampledir_rel+"/diagnostic/"+diag_out_file;
		}
		else{  // Embeds output into file
			diag_out_file = "[["+endli+di.te+"]]";
		}
		
		stringstream ss;
		ss << "diagnostics-inf chain=" << di.ch << " file=\""+diag_out_file+"\"" << endl;
		ss << endl;
		
		if(com_op == true) cout << ss.str();
		else fout << ss.str();
	}
}

			
/// Initialises structure for trans_diag observations
vector <TransDiagSpecies> Output::trans_diag_init() const
{
	auto T = model.details.T;

	vector <TransDiagSpecies> trans_diag;
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		TransDiagSpecies tdsp;
		tdsp.n = 0;
		tdsp.exp_num.resize(sp.tra_gl.size());
		for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
			tdsp.exp_num[tr].resize(T,0);
		}
		
		trans_diag.push_back(tdsp);
	}
	
	return trans_diag;
}


/// Adds results from particle onto h_bin
void Output::trans_diag_add(vector <TransDiagSpecies> &trans_diag, const vector <Particle> &part) const
{
	auto T = model.details.T;

	for(const auto &pa : part){
		for(auto p = 0u; p < model.nspecies; p++){
			const auto &sp = model.species[p];
		
			const auto &exp_num = pa.species[p].exp_num; 
			auto &exp_num_add = trans_diag[p].exp_num;
					
			trans_diag[p].n++;
			for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
				for(auto ti = 0u; ti < T; ti++){
					exp_num_add[tr][ti] += exp_num[tr][ti];
				}
			}
		}
	}
}	


/// Gets statistics from a vector
Stat Output::get_statistic(vector <double> &vec) const 
{
	auto n = vec.size();
	if(n == 0) emsg("Zero vector size");
	
	Stat stat;
	
	auto sum = 0.0, sum2 = 0.0;
	for(auto i = 0u; i < vec.size(); i++) {
		sum += vec[i];
		sum2 += vec[i]*vec[i];
	}
	sum /= n;
	sum2 /= n;

	stat.mean = sum;
	auto var = sum2 - sum*sum; if(var < 0) var = 0;
	stat.sd = sqrt(var);
	
	sort(vec.begin(),vec.end());

	if(n >= 2) {
		auto i = floor_int((n-1)*0.025);
		auto f = (n-1)*0.025 - i;
		stat.CImin = vec[i]*(1-f) + vec[i+1]*f;

		i = floor_int((n-1)*0.975);
		f = (n-1)*0.975 - i;
		stat.CImax = vec[i]*(1-f) + vec[i+1]*f;
	} 
	else{
		stat.CImin = vec[0];
		stat.CImax = vec[0];
	}

	return stat;
}


/// Takes an average for each generation
string Output::generation_average(string head, string content) const
{
	auto head_col = split(head,',');
	auto ncol = head_col.size()-1;
	
	Hash hash_col;
	
	vector <string> g_st;
	vector < vector < vector <double> > > vec;
	auto lines = split(content,'\n');
	for(auto i = 0u; i < lines.size(); i++){
		auto col_val = split(lines[i],',');
		
		if(col_val.size() == ncol+1){
			auto g = col_val[0];
			auto j = hash_col.find(g);
			if(j == UNSET){
				j = vec.size();
				hash_col.add(j,g);
				vector < vector <double> > vec_add;
				vec_add.resize(ncol);
				vec.push_back(vec_add);
				g_st.push_back(g);
			}
		
			for(auto k = 0u; k < ncol; k++){
				auto va = number(col_val[1+k]);
				vec[j][k].push_back(va);
			}
		}
	}		
	
	stringstream ss;
	ss << head;
	for(auto j = 0u; j < vec.size(); j++){
		ss << g_st[j];
		for(auto k = 0u; k < ncol; k++){
			auto stat = get_statistic(vec[j][k]);
			ss << "," << stat.mean << "|" << stat.CImin << "|" << stat.CImax;
		}	
		ss << endl;
	}
	return ss.str();
}


/// Numbers particles (if they are unset)
void Output::number_part(vector <Particle> &part) const
{
	auto num = 1;
	for(auto &pa : part){
		if(pa.s == UNSET){ pa.s = num; num++;}
	}
}


/// Outputs the final cpu time
void Output::final_time(long sec, long op_sec) const 
{
	auto op_st = " ("+tstr((int)((op_sec*100)/(sec+TINY)))+"% outputting)";
	
	cout << endl << "Total CPU time: ";
	
	if(sec < 60){
		cout << sec << " seconds" << op_st << endl;
		return;
	}
	
	cout << std::fixed;
  cout << setprecision(1);
		
	auto min = sec/60.0;
	if(min < 60){
		cout << min << " minutes" << op_st << endl;
		return;
	}
	
	auto hour = min/60.0;
	if(hour < 24){
		cout << hour << " hours" << op_st << endl;
		return;
	}
	
	auto day = hour/24.0;
	cout << day << " days" << op_st << endl;
}


/// Outputs the memory usage
void Output::final_memory_usage() const 
{
	auto mem = memory_usage();

#ifdef USE_MPI	
	mpi.sum(mem);
#endif
	
	if(op()){
		auto per = (unsigned int)(100*mem/total_memory());
		cout << "Total memory consumption: " << mem_print(mem);
		cout << " (" << per << "% of computer)" << endl;
	}
}
