// This file contains output functions

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
 
using namespace std;

#include "output.hh"
#include "utils.hh"

/// Initialises the output 
Output::Output(unsigned int _chain, const Model &model, const Input &input, Mpi &mpi) : model(model), mpi(mpi)
{
	chain = _chain;
	
	timer.resize(OUTTIMER_MAX); for(auto &ti : timer) ti = 0;
		
	outputdir = input.outputdir;
	
	ensure_directory(outputdir);                      // Creates the output directory
	
	if(op()){
		lines_raw = input.lines_raw;
		
		auto i = 0u;
		while(i < lines_raw.size())
		{
			auto st = trim(lines_raw[i]);
			
			auto spl = split(st,' ');
			auto command = spl[0];
			
			auto flag = false;
			
			if(st == "# OUTPUT") flag = true;
			
			if(command == "sim-param" || command == "sim-state"){
				flag = true;
			}
			
			if(command == "inf-param" || command == "inf-state"){
				if(model.mode != PPC) flag = true;
			}
			
			if(command == "post-sim-param" || command == "post-sim-state"){
				flag = true;
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
		
		lines_raw.push_back("");
		lines_raw.push_back("# OUTPUT");
		lines_raw.push_back("");
	}
	
	trace_init();
};


/// Create a directory if it doesn't already exist
void Output::ensure_directory(const string &path) const 
{
	struct stat st;
	if (stat(path.c_str(), &st) == -1){  	            // Directory not found

#ifdef WINDOWS
		int ret = mkdir(path.c_str());
#else	
		int ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
		
		if(ret == -1){
			emsg("Error creating the directory '"+path+"'");
		}
	}
}


/// Checks that the file has been successfully open
void Output::check_open(ofstream &fout, string file) const 
{
	if(!fout) emsg("Could not write to file '"+file+"'");
}


/// Outputs a summary of the model and data to a file
void Output::summary(const Model &model) const
{
	auto file = outputdir+"/model.txt";
	ofstream fout(file);
	check_open(fout,file);

	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		fout << "SPECIES: " << sp.name << "   type: ";
		switch(sp.type){
		case POPULATION: fout << "population-based"; break;
		case INDIVIDUAL: fout << "inidividual-based"; break;
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
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i < N; i++){
					auto &par = model.param[ieg.omega[j][i]];
					fout <<par.name << "  ";
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
					fout << "t=" << ob.t << " ";
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
		fout << "  ";
		fout << "auto: "; if(par.auto_value == true) fout << "true"; else fout << "false";
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
			for(auto val : par.spline_info.knot_times) fout << val << ",";
			fout << endl;
		}

		switch(model.mode){
		case SIM:
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
	auto file = outputdir+"/proposal.txt";
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
		
	case SPLINE_AFFECT:
		{
			ss << "AFFECT Spline affect " << model.spline[al.num].name << endl;	
		}
		break;
		
	case DIV_VALUE_AFFECT:
		{
			auto eq = model.species[al.num].markov_eqn[al.num2].eqn_ref;
			ss << "AFFECT Div Value " << model.eqn[eq].te_raw << endl;
		}
		break;
		
	case DIV_VALUE_FAST_AFFECT:
		{
			ss << "AFFECT Div Value Fast" << endl;
			auto eq = model.species[al.num].markov_eqn[al.num2].eqn_ref;
			ss << "AFFECT Div Value Fast" << model.eqn[eq].te_raw << endl;
		}
		break;
		
	case DIV_VALUE_LINEAR_AFFECT:
		{
			ss << "AFFECT Div Value Linear";
			{
				string str = "";
				for(auto k : al.linear_prop.me){
					auto eq = model.species[al.num].markov_eqn[k].eqn_ref;
					str += model.eqn[eq].te_raw +", ";
				}
				ss << trunc(str);
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
		ss << "AFFECT popnum ind w" << endl;
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
			else{
				ss << al.list.size() << " " << al.map.size() << ": ";
				for(auto t : al.list) ss << t << ",";

			}				
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
	case EXP_RATE: return "Exponential"; 
	case EXP_RATE_NM: return "Exponential NM"; 
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
	auto file = outputdir+"/data.txt";
	ofstream fout(file);
	check_open(fout,file);
	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		
		fout << "SPECIES: " << sp.name << endl;
		
		fout << "ENTER:" << endl;
		for(const auto &en : sp.enter){
			fout << en.name << " "<< en.time;
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
				fout << "," << ev.t << "  ";
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
			fout << "t=" << pd.t << "  value=" << pd.value;
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
			fout << "start=" << ptd.tmin << "  end=" << ptd.tmax;
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


/// Prints a set of parameters
void Output::print_param(const vector <double> &vec) const
{
	if(com_op == true) return;
	
	cout << "Parameters:" << endl;
	for(auto th = 0u; th < model.param_vec.size(); th++){
		cout << model.param_vec[th].name << " "<< vec[th] << endl;
	}
	cout << endl;
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
void Output::trace_init()
{
	param_out << "State"; 
	for(auto th = 0u; th < model.param_vec.size(); th++){
		const auto &par = model.param_vec[th];	
		
		if(model.param[par.th].trace_output == true){
			param_out << ",\"" << replace(par.name,"→","->") << "\"";
		}
	}

	for(auto i = 0u; i < model.derive.size(); i++){
		const auto &der = model.derive[i];
		if(der.time_dep == false){
			if(der.dep.size() == 0){
				param_out << ",\"" << der.name << "\"";
			}
			else{
				for(auto j = 0u; j < der.eq.size(); j++){
					param_out << ",\"" << der.name << "_";
					for(auto k = 0u; k < der.dep.size(); k++){
						const auto &dp = der.dep[k];
						auto m = (unsigned int)(j/dp.mult)%dp.list.size();
						if(k != 0) param_out << ",";
						param_out << dp.list[m];
					}
					param_out << "\"";
				}
			}
		}
	}
	
	if(model.trans_tree){
		param_out << ",N^origin,N^infected,N^mut-tree,N^mut-origin,N^unobs,t^root";
	}
	
	for(const auto &sp : model.species){
		if(sp.type == INDIVIDUAL) param_out << ",N^" << sp.name << "-total";;
	}
	
	ic_output_head();
	
	param_out << ",L^markov,L^non-markov,L^ie,L^dist,L^obs,L^genetic-proc,L^genetic-obs,L^init,Prior";

	param_out << endl;
}


/// Outputs the parameter file columns for initial conditions
void Output::ic_output_head() 
{
	for(const auto &sp : model.species){
		const auto &ic = sp.init_cond;
		if(ic.type == INIT_POP_DIST){
			if(ic.focal_cl == UNSET){
				param_out << ",N^" << sp.name;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					const auto &co = sp.comp_gl[c];
					if(co.erlang_hidden == false) param_out << ",f^init_(" << co.name << ")";
				}
			}
			else{
				const auto &claa = sp.cla[ic.focal_cl];
				for(auto c = 0u; c < claa.comp.size(); c++){
					const auto co = claa.comp[c];
					if(!co.erlang_hidden) param_out << ",N^init_(" << co.name << ")";
				}
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(cl != ic.focal_cl){
						const auto &claa = sp.cla[cl];
						for(auto c = 0u; c < claa.comp.size(); c++){
							const auto co = claa.comp[c];
							if(!co.erlang_hidden) param_out << ",f^init_(" << co.name << ")";
						}
					}
				}
			}
		}
	}
}
	
	
/// Outputs the initial conditions
void Output::ic_output(const Particle &part) 
{
	for(auto p = 0u; p < model.species.size(); p++){
		const auto &sp = model.species[p];
		const auto &ic = sp.init_cond;
		const auto &icv = part.species[p].init_cond_val;
		
		if(ic.type == INIT_POP_DIST){
			if(ic.focal_cl == UNSET){
				param_out << "," << icv.N_total;
				for(auto c = 0u; c < sp.comp_gl.size(); c++){
					const auto &co = sp.comp_gl[c];
					if(!co.erlang_hidden) param_out << "," << icv.frac[c];
				}
			}
			else{
				const auto &claa = sp.cla[ic.focal_cl];
				for(auto c = 0u; c < claa.comp.size(); c++){
					const auto &co = claa.comp[c];
					if(!co.erlang_hidden) param_out << "," << icv.N_focal[c];
				}
				for(auto cl = 0u; cl < sp.ncla; cl++){
					if(cl != ic.focal_cl){
						const auto &claa = sp.cla[cl];
						for(auto c = 0u; c < claa.comp.size(); c++){
							const auto &co = claa.comp[c];
							if(!co.erlang_hidden) param_out << "," << icv.frac_focal[cl][c];
						}
					}
				}
			}
		}
	}
}


/// Outputs the burnin fraction (needed for anneal scan)	
void Output::set_output_burnin(double burnin_frac)
{
	if(!op()) return;
	
	for(auto j = 0u; j < lines_raw.size(); j++){
		auto st = lines_raw[j];
		auto spl = split(st,' ');
		if(spl[0] == "inference"){
			auto k = 0u; while(k < st.length()-12 && st.substr(k,12) != " burnin-frac") k++;
			if(k < st.length()-12) st = st.substr(0,k);
			st += " burnin-frac=" + to_str(burnin_frac);
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
void Output::param_sample(unsigned int s, const State &state)
{
	timer[PARAM_OUTPUT] -= clock();
	
	auto part = state.generate_particle();
	param_sample(s,part);
	
	timer[PARAM_OUTPUT] += clock();
}


/// Outputs a parameter sample 
void Output::param_sample(unsigned int s, const Particle &part)
{
	param_out << s;
	for(auto th = 0u; th < model.param_vec.size(); th++){
		const auto &par = model.param_vec[th];
		if(model.param[par.th].trace_output == true){
			param_out << "," << part.param_val[th];
		}
	}
		
	auto dir_out = part.dir_out;
	for(auto i = 0u; i < model.derive.size(); i++){
		const auto &der = model.derive[i];
		if(der.time_dep == false){
			const auto &out = dir_out[i];
			
			if(der.dep.size() == 0){
				param_out << "," << out.value_str[0];
			}
			else{
				for(auto j = 0u; j < der.eq.size(); j++){
					param_out << "," << out.value_str[j];
				}
			}
		}
	}

	if(model.trans_tree){ 
		const auto &tts = part.trans_tree_stats;
		param_out << "," << tts.N_origin << "," << tts.N_inf << "," << tts.N_mut_tree << "," << tts.N_mut_origin << "," << tts.N_unobs << "," << tts.t_root;
	}
	
	for(auto p = 0u; p < model.species.size(); p++){
		if(model.species[p].type == INDIVIDUAL){
			param_out << "," << part.species[p].individual.size();
		}
	}
	
	ic_output(part);
	
	const auto &like = part.like;
	
	if(std::isnan(like.init_cond)) emsg("pp");
	param_out << "," << like.markov << "," << like.nm_trans << "," << like.ie << "," << like.dist << "," << like.obs << "," << like.genetic_process << "," << like.genetic_obs << "," << like.init_cond << "," << like.prior+like.spline_prior+like.init_cond_prior;
	
	param_out << endl;
}


/// Outputs a state sample
void Output::state_sample(unsigned int s, const State &state)
{
	timer[STATE_OUTPUT] -= clock();
	state_sample(s,state.generate_particle());
	timer[STATE_OUTPUT] += clock();
}

	
/// Outputs a state sample
void Output::state_sample(unsigned int s, const Particle &part)
{
	auto value = param_value_from_vec(part.param_val); 
	
	stringstream ss;
	ss << "<<STATE " << s << ">>" << endl << endl;
	
	ss << "<PARAMETERS>" << endl;
	
	ss << output_param(value,part);

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
					ss << replace(sp.tra_gl[i].name,"→","->") << ":";
					
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
				
				/*
				auto flag = false;
				for(auto i = 0u; i < sp.tra_gl.size(); i++){
					if(flag == true) ss << ","; 
					flag = true;
					ss << replace(sp.tra_gl[i].name,"→","->");
				}
				ss << endl;
					
				for(auto ti = 0u; ti < T; ti++){
					auto flag = false;
				
					for(auto j = 0u; j < ssp.trans_num.size(); j++){
						if(flag == true) ss << ",";
						flag = true;
						ss << ssp.trans_num[j][ti];
					}
					ss << endl;
				}
				ss << endl;
				*/
			}
			break;
			
		case INDIVIDUAL:
			{
				ss << "<INDIVIDUALS>" << endl;
				ss << "name,source";
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
						auto c = UNSET;//ind.cinit;
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
								}
								break;
								
							case M_TRANS_EV: case NM_TRANS_EV: 
								{
									const auto &trg = sp.tra_gl[ev.tr_gl];
									if(trg.i == UNSET) ss <<  replace(trg.name,"→","->");
									else{
										const auto &tr = sp.cla[trg.cl].tra[trg.tr];
										ss << replace(tr.name,"→","->");
									}
									
									const auto &iif = ev.ind_inf_from;
									if(iif.p != UNSET){
										ss << "[";
										if(iif.p == OUTSIDE_INF) ss << "OUT_INF"; 
										else{
											auto pp = iif.p;
											if(pp != p) ss << model.species[pp].name << "|";
											ss << part.species[pp].individual[iif.i].name;
										}
										ss << "]";
									}
									
									if(c != trg.i) emsg("Inconsist");
									cnew = trg.f;
								}
								break;
							}
							ss << ":" << ev.t;
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
		ss << "<PHYLOTREE>" << endl;
		//if(part.inf_node.size() == 0) emsg(" zero");
		for(const auto &in : part.inf_node){
			ss << part.species[in.p].individual[in.i].name  << "|" << in.t_start << "|";	
			
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
				
				ss << "|" << iev.index << "|" << iev.t << "|";
				if(iev.mut_num == UNSET) ss << "NA";
				else ss << iev.mut_num;
			}
			ss << endl;
		}
		ss << endl;
	}
	
	state_out << ss.str();
}


/// Generates the final output files
void Output::generate_files()
{
	string state_out_str = "";
	if(use_ind_key){
		stringstream ss;
		ss << "{" << endl;
		ss << "  timepoint " << model.details.t_start << ":" <<  model.details.dt << ":" << model.details.t_end << endl;
		for(auto i = 0u; i < ind_key.size(); i++){
			ss << "  " << i << ":" << ind_key[i] << endl;
		}
		ss << "}" << endl << endl;
		state_out_str += ss.str();
	}
	state_out_str += state_out.str();
	
	/*
	if(com_op == true){
		cout << "<<PARAMS>>" << endl;
		cout << param_out.str();
		cout << "<<END>>" << endl;
		
		cout << "<<STATES>>" << endl;
		cout << state_out_str;
		cout << "<<END>>" << endl;
	}
	*/
	
	auto dir = outputdir+"/Sample";

	ensure_directory(dir);

	string ch = "1"; if(chain != UNSET) ch = tstr(chain);
	
	string add; if(chain != UNSET) add = "_"+tstr(chain);
	
	auto param_out_file = "param"+add+".csv";

	auto state_out_file = "state"+add+".csv";
	
	ofstream pout(dir+"/"+param_out_file);
	check_open(pout,param_out_file);
	pout << param_out.str();
	
	ofstream sout(dir+"/"+state_out_file);
	check_open(sout,state_out_file);
	sout << state_out_str;
	
	param_out_file = "Sample/"+param_out_file;
	state_out_file = "Sample/"+state_out_file;
	
	if(true){ // Embeds output into file
		param_out_file = "[["+endli+param_out.str()+"]]";
		state_out_file = "[["+endli+state_out_str+"]]";
	}
	
	switch(model.mode){
	case SIM:
		lines_raw.push_back("sim-param file=\""+param_out_file+"\"");
		lines_raw.push_back("");
		lines_raw.push_back("sim-state file=\""+state_out_file+"\"");
		break;
	
	case INF:
		lines_raw.push_back("inf-param chain="+ch+" file=\""+param_out_file+"\"");
		lines_raw.push_back("");
		lines_raw.push_back("inf-state chain="+ch+" file=\""+state_out_file+"\"");
		break;
		
	case PPC:
		lines_raw.push_back("post-sim-param file=\""+param_out_file+"\"");
		lines_raw.push_back("");
		lines_raw.push_back("post-sim-state file=\""+state_out_file+"\"");
		break;
		
	default: emsg("Should not be default12"); break;
	}
}


/// Generates values from parameter vector
vector < vector <double> > Output::param_value_from_vec(const vector <double> &param_val) const 
{
	vector < vector <double> > value;
	
	value.resize(model.param.size());
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		if(par.variety != CONST_PARAM || par.factor){
			value[th].resize(par.N,UNSET);
		}
	}
	
	if(param_val.size() != model.param_vec.size()) emsg("param_val the wrong size");
	
	for(auto i = 0u; i < model.param_vec.size(); i++){
		const auto &mpv = model.param_vec[i];
		value[mpv.th][mpv.index] = param_val[i]; 
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
			for(auto i = 0u; i < par.N; i++) value[th][i] = par.value[i].value;
		}
	}

	return value;
}


/// Outputs a string version of parameters which are either constant or non-constant 
string Output::output_param(const vector < vector <double> > &value, const Particle &part) const
{
	stringstream ss;
	
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		
		if(par.variety != CONST_PARAM && par.trace_output == true){
			ss << "\"" << replace(par.name,"→","->") << "\"";
			if(par.dep.size() == 0){
				ss << "," << value[th][0] << endl;
			}
			else{
				ss << endl;
				ss << "  ";
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
	
	const auto &like = part.like;
	ss << "L^markov," << like.markov << endl;
	ss << "L^non-markov," << like.nm_trans << endl;
	ss << "L^ie," << like.ie << endl;
	ss << "L^dist," << like.dist << endl;
	ss << "L^obs," << like.obs << endl;
	ss << "L^init," << like.init_cond << endl;
	ss << "Prior," << like.prior << endl;
	
	ss << endl;
	
	return ss.str();
}


/// Outputs an updated version of the input file
void Output::updated_file(string file)
{
	if(com_op == true){
		cout << "<<OUTPUT FILE>>" << endl;
		if(chain == 0){
			for(auto li : lines_raw){
				cout << li << endl;
			}
		}
		else{ // For other chains just shows the output
			auto i = 0u; 
			while(i < lines_raw.size() && lines_raw[i] != "# OUTPUT") i++;
			while(i < lines_raw.size()){
				cout << lines_raw[i] << endl;
				i++;
			}
		}
		cout << "<<END>>" << endl;
		return;
	}
	
#ifdef USE_MPI	
	mpi.transfer_lines_raw(lines_raw);
#endif
	
	if(op()){
		ofstream fout(file);
	
		for(auto li : lines_raw){
			fout << li << endl;
		}
	}	
}


/// Prints out observation events 
string Output::print_obs_data(unsigned int p, const vector <ObsData> &obs) const 
{
	stringstream ss;
	
	const auto &sp = model.species[p];
	
	ss << "Observations:" << endl;
	for(auto &ob : obs){
		const auto &claa =  sp.cla[ob.cl];
		
		ss << ob.t << " ";
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
