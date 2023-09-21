// This file defines the data

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
Output::Output(unsigned int chain, const Model &model, const Input &input) : model(model)
{
	outputdir = input.outputdir;
	
	ensure_directory(outputdir);   // Creates the output directory
	lines_raw = input.lines_raw;
	
	auto i = 0u;
	while(i < lines_raw.size())
	{
		auto st = trim(lines_raw[i]);
		
		auto spl = split(st,' ');
		auto command = spl[0];
		
		if(st == "# OUTPUT" || command == "sim-param" || command == "sim-state" || command == "inf-param" || command == "inf-state"){
			//command == "sim-const" || command == "inf-const"
			lines_raw.erase(lines_raw.begin()+i); 
		}
		else i++;
	}
	
	while(lines_raw.size() > 0 && lines_raw[lines_raw.size()-1] == ""){
		lines_raw.erase(lines_raw.begin()+lines_raw.size()-1); 
	}
	
	lines_raw.push_back("");
	lines_raw.push_back("# OUTPUT");
	lines_raw.push_back("");
	
	trace_init(chain);
};


/// Create a directory if it doesn't already exist
void Output::ensure_directory(const string &path) const 
{
	//if(false) cout << path;
	
	struct stat st;
	if (stat(path.c_str(), &st) == -1){  	// Directory not found

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

		for(const auto &cl : sp.cla){
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
				if(tr.i == SOURCE) fout << "+"; else fout << cl.comp[tr.i].name;

				fout << " -> ";
				if(tr.f == SINK) fout << "-"; else fout << cl.comp[tr.f].name;

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
			/*
			fout << "Branch points:" << endl;
			for(const auto &bp : cl.branch_point){
				fout << "Compartment: " << cl.comp[bp.c].name << "   Trans: ";
				for(const auto tr : bp.tra) fout << cl.tra[tr].name << ", ";
				fout << endl;
			}
			*/
			fout << endl;
		}
	
		if(true){
			fout << "Global Compartments:" << endl;
			for(const auto &co_gl : sp.comp_gl){
				fout << co_gl.name << endl;
				for(auto cl = 0u; cl < sp.ncla; cl++){
					const auto &cla = sp.cla[cl];
					const auto &ti = co_gl.trainfo[cl];
					
					fout << cla.name << " ";
					if(ti.branch == true) fout << "branch"; else fout << "no branch";
					fout << ": ";
					for(auto tref : ti.tra_ref){
						fout << sp.tra_gl[tref.tr].name << ", ";
						auto th = tref.bp_th;
						if(th != UNSET) fout << "BP: " << model.param_vec[th].name << " ";
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
				if(tr_gl.i == SOURCE) fout << "+";
				else fout << sp.comp_gl[tr_gl.i].name;
				fout << " -> ";
				if(tr_gl.f == SINK) fout << "-";
				else fout << sp.comp_gl[tr_gl.f].name;
				fout << endl;
			}
		}

		fout << "DATA:" << endl;
		
		for(const auto &ds : sp.source){
			switch(ds.cname){
			case INIT_POP_SIM: fout << "Init. pop sim"; break;
			case ADD_IND_SIM: fout << "Add ind sim"; break;
			case REMOVE_IND_SIM: fout << "Remove ind sim"; break;
			case MOVE_IND_SIM: fout << "Move ind sim"; break; 
			case INIT_POP: fout << "Init pop"; break; 
			case ADD_IND: fout << "Add ind"; break; 
			case REMOVE_IND: fout << "Remove ind"; break;
			case MOVE_IND: fout << "Move ind"; break; 
			case INIT_POP_PRIOR: fout << "Init pop prior"; break; 
			case COMP_DATA: fout << "Comp data"; break; 
			case TRANS_DATA: fout << "Trans data"; break;
			case SOURCE_DATA: fout << "Source data"; break; 
			case SINK_DATA: fout << "Sink data0"; break;
			case TEST_DATA: fout << "Test data"; break; 
			case POP_DATA: fout << "Pop data"; break; 
			case POP_TRANS_DATA: fout << "Pop trans data"; break; 
			case SET_TRAPS_DATA: fout << "Set traps data"; break; 
			case IND_TRAPPED_DATA: fout << "Ind trapped data"; break;
			case GENETIC_DATA: fout << "Genetic data"; break;
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
			if(nmt.bp_eq_ref != UNSET) fout << model.eqn[nmt.bp_eq_ref].te_raw << "  ";
			for(auto eq : nmt.dist_param_eq_ref) fout << model.eqn[eq].te_raw << "  ";
			fout << endl;
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
	
	if(true){
		fout << "BRANCH PARAM GROUPS:" << endl;
		for(const auto &bg : model.branch_param){
			for(auto th : bg.group) fout << model.param_vec[th].name << ", ";
			fout << endl;				
		}
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
		
		fout << " ME: ";
		for(const auto &mer : po.markov_eqn_ref){
			const auto &sp = model.species[mer.p];
			auto eq = sp.markov_eqn[mer.e].eqn_ref;
			fout << model.eqn[eq].te_raw << "  ";
		}
		fout << endl;
	}
		
	fout << endl;
		/*
	
		fout << "Capture events:" << endl;
		for(const auto &ce : capev){
			fout << ce.name << " " << cla[ce.cl].name << "  ";
			switch(ce.type){
				case CAPEV_TRANS: fout << "TRANS " << cla[ce.cl].trans[ce.tr].name << " " << ce.tmin << " - " << ce.tmax; break;
				default: emsg("Model output error 3a"); break;
			}
			fout << endl;
		}		
		
		fout << endl;
		
		fout << "Capture event timeline:" << endl;
		for(const auto &cetl : capev_timeline){
			fout << cetl.time << "   ";
			switch(cetl.type){
				case CAPEV_START:
					fout << "START ";
					for(auto tr : cetl.tr_list) fout << trans[tr].name << ","; 
					break;
					
				case CAPEV_END:
					fout << "END ";
					for(auto tr : cetl.tr_list) fout << trans[tr].name << ","; 
					break;
				
				default: emsg("Capev problem"); break;
			}
			fout << endl;
		}
		
		fout << endl;
		
		fout << "Indivuduals:" << endl;
		for(const auto &ind : individual){
			fout << ind.id << ":    ";
			for(const auto &ob : ind.obs){
				fout << "t=" << ob.t << " ";
				switch(ob.type){
					case OBS_STATE: fout << cla[ob.cl].comp[ob.value].name; break;
					case OBS_TRANS: fout << cla[ob.cl].trans[ob.value].name; break;
					default: emsg("Model output error 4"); break;
				}
				fout << "    ";
			}
			fout << endl;
		}
		
		*/
	

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
		fout << "age_dep: "; if(par.age_dep == true) fout << "true"; else fout << "false";
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

		case INF:
			break;

		default:  emsg("Model output error10"); break;
		}
		
		fout << endl;
		fout << endl;
	}
	
	fout << "PARAMETER VECTOR:" << endl;
	for(const auto &par : model.param_vec){
		fout << par.name << "   affect: " << endl;
		for(auto &al : par.affect_like){
			fout << print_affect_like(al);
		}
		fout << endl;
	}
	fout << endl;
	
	if(true){
		fout << "SPLINES:" << endl;
		for(const auto &spl : model.spline){
			fout << spl.name << endl;
		}
		fout << endl;
	}
	
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
	
	if(true){
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


// Prints how the state properties are altered
string Output::print_affect_like(const AffectLike &al) const
{
	stringstream ss;
	
	switch(al.type){
	case LIKE_OBS_AFFECT:
		{
			ss << "AFFECT like obs for species " << model.species[al.num].name;
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
			const auto &nm = model.species[al.num].nm_trans[al.num2];
			ss << "AFFECT Nonmarkovian trans affect: " << model.eqn[nm.dist_param_eq_ref[0]].te_raw << endl;	
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
	}
	
	if(al.list.size() != 0){
		if(al.type == POP_AFFECT){
			ss << "Populations affected: ";	
			for(auto k : al.list) ss << model.pop[k].name << ","; 
			ss << endl;	
		}
		else{		
			ss << "Times affected: ";	
			if(al.list.size() == al.map.size()) ss << " All";
			else for(auto t : al.list) ss << t << ","; 
			ss << endl;	
		}
	}

	return ss.str();	
}


/// Converts a type into text;
string Output::transtype_text(TransType type) const 
{
	switch(type){
	case EXP_RATE: return "Exp(rate)"; 
	case EXP_MEAN: return "Exp(mean)";
	case GAMMA: return "Gamma";
	case ERLANG: return "Erlang"; 
	case LOG_NORMAL: return "Log-normal"; 
	case WEIBULL: return "Weibull";
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
		for(const auto &ind : sp.individual){
			fout << ind.name << ": ";
			for(const auto &ev : ind.ev){
				switch(ev.type){
				case ENTER_EV: fout << "enter " << sp.comp_gl[ev.c].name; break;
				case LEAVE_EV: fout << "leave"; break;
				case MOVE_EV: fout << "move to " << sp.cla[ev.cl].comp[ev.c].name; break;
				case M_TRANS_EV: fout << "markov trans " << sp.cla[ev.cl].tra[ev.tr].name; break;	
				case NM_TRANS_EV: fout << "non-markov trans " << sp.cla[ev.cl].tra[ev.tr].name; break;	
				}
				fout << "," << ev.t << "  ";
			}
			fout << endl;
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
		const auto &ic = initc[p];
		
		const auto &sp = model.species[p];
		
		cout << "Initial conditions for '"+sp.name+"':" << endl;
		
		switch(ic.type){
		case POP_INIT:
			for(auto c = 0u; c < sp.comp_gl.size(); c++){
				cout << sp.comp_gl[c].name << " " << ic.cpop[c] << endl;
			}
			break;
			
		case ZEROPOP_INIT: 
			cout << "zero pop\n";
			break;
			
		default: emsg("Should not be default11"); break;
		}
	}
	cout << endl;
}


/// Outputs a parameter sample 
void Output::trace_init(unsigned int chain)
{
	string add; if(chain != UNSET) add = "_"+tstr(chain);
	
	param_out_file = "param"+add+".csv";

	state_out_file = "state"+add+".csv";
		
	string ch = "1"; if(chain != UNSET) ch = tstr(chain);
	
	switch(model.mode){
	case SIM:
		lines_raw.push_back("sim-param file=\"Sample/"+param_out_file+"\"");
		lines_raw.push_back("sim-state file=\"Sample/"+state_out_file+"\"");
		break;
	
	case INF:
		lines_raw.push_back("inf-param chain="+ch+" file=\"Sample/"+param_out_file+"\"");
		lines_raw.push_back("inf-state chain="+ch+" file=\"Sample/"+state_out_file+"\"");
		break;
		
	default: emsg("Should not be default12"); break;
	}
	
	param_out << "State";
	for(auto th = 0u; th < model.param_vec.size(); th++){
		const auto &par = model.param_vec[th];
		param_out << ",\"" << replace(par.name,"→","->") << "\"";
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
	
	param_out << ",L^markov,L^non-markov,L^ie,L^dist,L^obs,Prior";

	param_out << endl;
}


/// Outputs a parameter sample 
void Output::param_sample(unsigned int s, const State &state)
{
	auto part = state.generate_particle();
	param_sample(s,part);
}


/// Outputs a parameter sample 
void Output::param_sample(unsigned int s, const Particle &part)
{
	param_out << s;
	for(auto th = 0u; th < model.param_vec.size(); th++){
		param_out << "," << part.param_val[th];
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


	const auto &like = part.like;
	
	param_out << "," << like.markov << "," << like.nm_trans << "," << like.ie << "," << like.dist << "," << like.obs << "," << like.prior+like.spline_prior;

	param_out << endl;
}


/// Outputs a state sample
void Output::state_sample(unsigned int s, const State &state)
{
	state_sample(s,state.generate_particle());
}

	
/// Outputs a state sample
void Output::state_sample(unsigned int s, const Particle &part)
{
	auto value = param_value_from_vec(part.param_val); 
	
	stringstream ss;
	ss << "<<STATE " << s << ">>" << endl << endl;
	
	ss << "<PARAMETERS>" << endl;
	
	ss << output_param(value,part);

	ss << "<TIMEPOINT>" << endl;
	const auto &tp = model.timepoint;
	for(auto ti = 0u; ti < tp.size(); ti++){
		if(ti != 0) ss << ",";
		ss << tp[ti];
	}
	ss << endl << endl;
	
	for(auto p = 0u; p < model.nspecies; p++){
		const auto &sp = model.species[p];
		const auto &ssp = part.species[p];
	
		ss << "<SPECIES \"" << sp.name << "\">" << endl;
		ss << endl;
			
		if(model.species[p].type == POPULATION){
			ss << "<INITIAL POPULATION>" << endl;
			ss << "compartment,population" << endl;
			for(auto c = 0u; c < sp.comp_gl.size(); c++){
				if(sp.comp_gl[c].erlang_hidden == false){
					ss << sp.comp_gl[c].name << "," << ssp.cpop_init[c] << endl;
				}
			}
			ss << endl;
			
			ss << "<TRANSITIONS>" << endl;
			auto flag = false;
			for(auto i = 0u; i < sp.tra_gl.size(); i++){
				if(sp.tra_gl[i].erlang_hidden == false){
					if(flag == true) ss << ","; 
					flag = true;
					ss << replace(sp.tra_gl[i].name,"→","->");
				}
			}
			ss << endl;
			
			const auto &tp = model.timepoint;
			
			for(auto ti = 0u; ti < tp.size()-1; ti++){
				auto flag = false;
			
				for(auto j = 0u; j < ssp.trans_num.size(); j++){
					if(sp.tra_gl[j].erlang_hidden == false){
						if(flag == true) ss << ",";
						flag = true;
						ss << ssp.trans_num[j][ti];
					}
				}
				ss << endl;
			}
			ss << endl;
		}
		else{
			ss << "<INDIVIDUALS>" << endl;
			ss << "name,init";
			for(auto i = 0u; i < sp.ind_effect.size(); i++) ss << "," << sp.ind_effect[i].name;
			ss << ",events" << endl;
			for(auto i = 0u; i < ssp.individual.size(); i++){
				const auto &ind = ssp.individual[i];
				ss << ind.name << ",";
				if(ind.cinit == SOURCE) ss << "+";
				else{
					if(ind.cinit != UNSET) ss << sp.comp_gl[ind.cinit].name;
				}
				
				for(auto i = 0u; i < sp.ind_effect.size(); i++) ss << "," << ind.exp_ie[i]; 
				ss << ",";
				auto c = ind.cinit;
				for(auto e = 0u; e < ind.ev.size(); e++){
					const auto &ev = ind.ev[e];
					auto cnew = c;
					
					if(e != 0) ss << " ";
					switch(ev.type){
					case ENTER_EV:
						if(c != UNSET) emsg("Inconsistent");
						cnew = ev.c_gl; 
						ss << sp.comp_gl[cnew].name;
						break;

					case LEAVE_EV: 
						cnew = UNSET; 
						ss << "-"; 
						break;
						
					case MOVE_EV: 
						if(c != UNSET){
							auto cl = ev.tr_gl;
							cnew = sp.update_c_comp(c,cl,ev.c_gl);
							ss << "move->" << sp.cla[cl].comp[ev.c_gl].name;
						}
						break;
						
					case M_TRANS_EV: case NM_TRANS_EV: 
						const auto &trg = sp.tra_gl[ev.tr_gl];
						if(trg.i == SOURCE) ss <<  replace(trg.name,"→","->");
						else{
							const auto &tr = sp.cla[trg.cl].tra[trg.tr];
							ss << replace(tr.name,"→","->");
						}
						
						if(c != trg.i) emsg("Inconsist");
						cnew = trg.f;
						break;
					}
					ss << ":" << ev.t;
					c = cnew;
				}
				ss << endl;
			}
			ss << endl;
		}
	}
	
	if(model.derive.size() > 0){
		auto dir_out = part.dir_out;//state.derive_calculate();
		
		ss << "<DERIVED>" << endl;
		for(auto i = 0u; i < model.derive.size(); i++){
			const auto &der = model.derive[i];
			const auto &out = dir_out[i];
			
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
	
	state_out << ss.str();
}


/// Generates the final output files
void Output::generate_files()
{
	if(com_op == true){
		cout << "<<PARAMS>>" << endl;
		cout << param_out.str();
		cout << "<<END>>" << endl;
		
		cout << "<<STATES>>" << endl;
		cout << state_out.str();
		cout << "<<END>>" << endl;
	}
	
	auto dir = outputdir+"/Sample";

	ensure_directory(dir);

	ofstream pout(dir+"/"+param_out_file);
	check_open(pout,param_out_file);
	pout << param_out.str();
	
	ofstream sout(dir+"/"+state_out_file);
	check_open(sout,state_out_file);
	sout << state_out.str();
}


/// Generates parameters values from parameter vector
vector < vector <double> > Output::param_value_from_vec(const vector <double> &param_val) const 
{
	vector < vector <double> > value;
	
	value.resize(model.param.size());
	for(auto th = 0u; th < model.param.size(); th++){
		const auto &par = model.param[th];
		if(par.variety != CONST_PARAM){
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
		
		if(par.variety != CONST_PARAM){
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
	ss << "Prior," << like.prior << endl;
	
	ss << endl;
	
	return ss.str();
}
	

/// Outputs an updated version of the input file
void Output::updated_file(string file) const 
{
	if(com_op == true) return;
		
	ofstream fout(file);
	
	for(auto li : lines_raw){
		fout << li << endl;
	}
}


