// Initialises data structure from sources

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
 
using namespace std;

#include "utils.hh"
#include "species.hh"

/// Initialises the data sources
void Species::initialise_data()
{
	T = timepoint.size()-1;
	
	trans_not_allow.resize(T); pop_trans_ref.resize(T);
	for(auto ti = 0u; ti < T; ti++){
		trans_not_allow[ti].resize(tra_gl.size(),false);
		pop_trans_ref[ti].resize(tra_gl.size());
	}
	
	init_cond.type = UNSET_INIT;

	// Adds a transition filter on to transition data source
	for(auto &so : source){
		switch(so.cname){
		case TRANS_DATA: case SOURCE_DATA: case SINK_DATA:
			{	
				vector <bool> filt;
				for(auto tr = 0u; tr < tra_gl.size(); tr++){
					const auto &trg = tra_gl[tr];
					if(trg.cl == so.cl && trg.tr == so.tr) filt.push_back(true);
					else filt.push_back(false);
				}
				so.trans_filt = filt;
			}
			break;
		
		default: break;
		}
	}
	
	for(const auto &so : source){
		switch(so.cname){
		case INIT_POP: init_pop_data(so); break;
		case ADD_IND: add_ind_data(so); break;
		case REMOVE_IND: remove_ind_data(so); break;
		case MOVE_IND: move_ind_data(so); break;
		case TRANS_DATA: case SOURCE_DATA: case SINK_DATA: trans_data(so); break;
		case COMP_DATA: comp_data(so); break;
		case TEST_DATA: test_data(so); break;
		case POP_DATA: population_data(so); break;
		case POP_TRANS_DATA: population_trans_data(so); break;
		default:
			emsg("data type not added:"); break;
		}
	}
	
	X_vector_order();
	
	if(init_cond.type == UNSET_INIT){
		init_cond.type = ZEROPOP_INIT;
		init_cond.cpop.resize(comp_gl.size(),0);
	}		
	
	order_data_events();   
}


/// init-pop command
void Species::init_pop_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	if(init_cond.type != UNSET_INIT) alert_source("Cannot set initial conditions twice",so);
	
	init_cond.type = POP_INIT;
	const auto &comp = comp_gl;
	
	auto foc_cl = so.focal_cl; 
	if(foc_cl != UNSET){
		vector < vector <double> > per;
		per.resize(ncla);
		for(auto cl = 0u; cl < ncla; cl++){
			for(auto c = 0u; c < cla[cl].ncomp; c++){
				double val = UNSET; if(cla[cl].comp[c].erlang_hidden == true) val = 0;
				per[cl].push_back(val);
			}
		}					
		
		if(tab.ncol != 2) alert_source("Should have two columns",so);
		for(auto r = 0u; r < tab.nrow; r++){
			auto name = tab.ele[r][0];
			
			auto flag = false;
			for(auto cl = 0u; cl < ncla; cl++){
				for(auto c = 0u; c < cla[cl].ncomp; c++){
					if(cla[cl].comp[c].name == name){
						if(per[cl][c] != UNSET){
							alert_source("The comartment '"+name+"' is set more than once",so,0,r);
						}					
						
						auto val = tab.ele[r][1];
						if(cl == foc_cl){
							per[cl][c] = number(val);
							if(per[cl][c] == UNSET) alert_source("The value '"+val+"' is not a number",so,1,r);
						}
						else{
							if(is_percent(val) == false){
								alert_source("The value '"+val+"' is not a percentage",so,1,r);
							}
							else{
								per[cl][c] = number(val.substr(0,val.length()-1));
							}
						}
						
						flag = true;
						break;
					}
				}
				if(flag == true) break;
			}
			
			if(flag == false){
				alert_source("The compartment '"+name+"' is not recognised",so,0,r);
			}
		}		

		// Fills in missing percentage
		for(auto cl = 0u; cl < ncla; cl++){
			const auto &claa = cla[cl];
			auto sum = 0.0;
			vector <unsigned int> list;
			
			for(auto c = 0u; c < claa.ncomp; c++){
				if(per[cl][c] == UNSET) list.push_back(c);
				else sum += per[cl][c];
			}
			
			string st = ""; 
			for(auto j = 0u; j < list.size(); j++){
				if(j != 0) st += ", ";
				st += claa.comp[list[j]].name;
			}
			
			if(cl == foc_cl){
				if(list.size() != 0){
					alert_source("Populations for the compartment(s) '"+st+"' must be set",so);
				}
			}
			else{
				if(list.size() == 0){
					alert_source("Not all percentages for compartmental population percentages in '"+claa.name+"' should be set",so);
				}
				
				if(list.size() > 2){
					alert_source("All but one of the compartmental population percentages in '"+st+"' should be set",so);
				}
				
				if(list.size() == 1){
					if(sum > 100.0){
						alert_source("Population percentages in '"+claa.name+"' add to over 100%",so);
					}
					else{
						per[cl][list[0]] = 100-sum; 
					}
				}
			}
		}	

		if(false){
			for(auto cl = 0u; cl < ncla; cl++){
				const auto &claa = cla[cl];	
				for(auto c = 0u; c < claa.ncomp; c++){
					cout << claa.comp[c].name << " " << per[cl][c] << " Percent\n";
				}
			}
		}
		
		init_cond.cpop.resize(comp.size(),0);
		for(auto c = 0u; c < comp.size(); c++){
			const auto &co = comp[c];
			
			auto fac = 1.0;
			for(auto cl = 0u; cl < ncla; cl++){
				if(cl == foc_cl){
					fac *= per[cl][co.cla_comp[cl]];
				}
				else{
					fac *= per[cl][co.cla_comp[cl]]/100;
				}
			}
			init_cond.cpop[c] = fac;
		}
	}
	else{
		if(tab.ncol != ncla+1){
			alert_source("Does not have the right number of columns",so);
			return;
		}
		
		init_cond.cpop.resize(comp.size(),UNSET);
		for(auto r = 0u; r < tab.nrow; r++){
			auto c = 0u;
			for(auto cl = 0u; cl < ncla; cl++){
				auto name = tab.ele[r][cl];
			
				const auto &claa = cla[cl];
				auto j = 0u; while(j < claa.ncomp && claa.comp[j].name != name) j++;
				if(j == claa.ncomp){ alert_source("Could not find '"+name+"'",so); return;}
				
				c += comp_mult[cl]*j;
			}
			
			auto ele = tab.ele[r][ncla];
			auto num = number(ele);
			if(num == UNSET){ alert_source("Not a number '"+ele+"'",so); return;}
			init_cond.cpop[c] = num;
		}
		
		for(auto i = 0u; i < comp.size(); i++){
			if(init_cond.cpop[i] == UNSET){
				alert_source("Population for '"+comp[i].name+"' not set",so); return;
			}
		}
	}
	
	if(false){
		cout << "init pop\n";
		for(auto c = 0u; c < comp.size(); c++){
			const auto &co = comp[c];
			cout << co.name << " " << init_cond.cpop[c] << " pop" << endl;
		}
	}
}
		
		
/// add-ind command
void Species::add_ind_data(const DataSource &so)
{
	const auto &tab = so.table;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		
		auto t_str = tab.ele[j][1];
		double t = number(t_str);
		if(t_str == "start") t = details.t_start;
		
		if(ncla+2 != tab.ncol) emsg("Columns not right");
		
		auto c_gl = 0u;
		for(auto cl = 0u; cl < ncla; cl++){
			auto val = tab.ele[j][cl+2];
			auto c = find_c(cl,val);
			if(c == UNSET) alert_source("Value '"+val+"' is not a compartment",so,cl+2,j);
			
			c_gl += c*comp_mult[cl];
		}
		
		EventData ev; 
		ev.type = ENTER_EV;
		ev.c = c_gl;
		ev.cl = UNSET;
		ev.tr = UNSET;
		ev.t = t;
		ind.ev.push_back(ev);
	}		
}


/// remove-ind command
void Species::remove_ind_data(const DataSource &so)
{
	const auto &tab = so.table;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto t_str = tab.ele[j][1];
		double t = number(t_str);
		if(t_str == "start") t = details.t_start;
		
		if(tab.ncol != 2) emsg("Columns not right");
		
		EventData ev; 
		ev.type = LEAVE_EV;
		ev.c = UNSET;
		ev.cl = UNSET;
		ev.tr = UNSET;
		ev.t = t;
		ind.ev.push_back(ev);
	}		
}


/// move-ind command
void Species::move_ind_data(const DataSource &so)
{
	const auto &tab = so.table;
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto t_str = tab.ele[j][1];
		double t = number(t_str);
		if(t_str == "start") t = details.t_start;
		
		if(tab.ncol != 3) emsg("Columns not right");
		
		auto cl = so.cl;
		auto val = tab.ele[j][2];
		auto c = find_c(cl,val);
		if(c == UNSET) alert_source("Value '"+val+"' is not a compartment",so,2,j);
			
		EventData ev; 
		ev.type = MOVE_EV;
		ev.c = c;
		ev.cl = cl;
		ev.tr = UNSET;
		ev.t = t;
		ind.ev.push_back(ev);
	}		
}


/// trans-data command
void Species::trans_data(const DataSource &so)
{
	const auto &tab = so.table;

	// Adds the times of events
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto val = tab.ele[j][1];
		if(val != "no"){
			auto t = number(val);
			if(t == UNSET){
				alert_source("Value '"+val+"' is not a number",so,1,j);
				return;
			}
			
			ObsData ob; 
			ob.so = so.index;
			ob.type = OBS_TRANS_EV;
			ob.c = UNSET;
			ob.cl = so.cl;
			ob.tr = so.tr;
			ob.t = t;
			ind.obs.push_back(ob);
		}
	}		
	
	
	// Adds the times of non-events
	for(auto j = 0u; j < tab.nrow; j++){
		//auto i = find_individual(tab.ele[j][0]);
		//auto &ind = individual[i];
			
		auto ti_min = 0u, ti_max = T; 
		
		switch(so.time_range){
		case ALL_TIME:
			{
				if(tab.ncol != 2) emsg("Columns not right");
			}
			break;
			
		case SPEC_TIME: 
			{
				if(tab.ncol != 2) emsg("Columns not right");
				ti_min = get_ti(so.time_start);
				ti_max = get_ti(so.time_end);
			}
			break;
	
		case FILE_TIME:
			{
				if(tab.ncol != 4) emsg("Columns not right");
				auto val = tab.ele[j][2];
				auto tmin = number(val);
				if(tmin == UNSET){
					alert_source("Value '"+val+"' is not a number",so,2,j);
					return;
				}
			
				val = tab.ele[j][3];
				auto tmax = number(val);
				if(tmax == UNSET){
					alert_source("Value '"+val+"' is not a number",so,3,j);
					return;
				}
				
				ti_min = get_ti(tmin);
				ti_max = get_ti(tmax);
			}
			break;
		}
		
		for(auto ti = ti_min; ti < ti_max; ti++){
			for(auto tr = 0u; tr < tra_gl.size(); tr++){
				if(so.trans_filt[tr] == true) trans_not_allow[ti][tr] = true;
			}
		}
	}
}


/// comp-data command
void Species::comp_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto val = tab.ele[j][1];
		auto t = number(val);
		if(t == UNSET){
			alert_source("Value '"+val+"' is not a number",so,1,j);
			return;
		}
		
		auto cl = so.cl;
		
		auto c_st = tab.ele[j][2];
		auto c = find_c(cl,c_st);
		if(c == UNSET){
			alert_source("Value '"+val+"' is not a compartment",so,2,j);
			return;
		}
		
		ObsData ob; 
		ob.so = so.index;
		ob.type = OBS_COMP_EV;
		ob.c = c;
		ob.cl = so.cl;
		ob.tr = UNSET;
		ob.t = t;
		ind.obs.push_back(ob);
	}		
}


/// test-data command
void Species::test_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	for(auto j = 0u; j < tab.nrow; j++){
		auto i = find_individual(tab.ele[j][0]);
		auto &ind = individual[i];

		auto val = tab.ele[j][1];
		auto t = number(val);
		if(t == UNSET){
			alert_source("Value '"+val+"' is not a number",so,1,j);
			return;
		}
		
		ObsData ob; 
		ob.type = OBS_TEST_EV;
		ob.so = so.index;
		ob.c = UNSET;
		ob.cl = so.cl;
		ob.tr = UNSET;
		ob.t = t;
		ob.Se_obs_eqn = add_to_vec(obs_eqn,so.obs_model.Se.eq_ref);
		ob.Sp_obs_eqn = add_to_vec(obs_eqn,so.obs_model.Sp.eq_ref);
		
		auto res = tab.ele[j][2];
		if(res == so.obs_model.diag_pos) ob.test_res = true;
		else{
			if(res == so.obs_model.diag_neg) ob.test_res = false;
			else{
				alert_source("Value '"+res+"' is not a positive or negative test result",so,2,j);
				return;
			}
		}
		
		ind.obs.push_back(ob);
	}		
}


/// pop-data command
void Species::population_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	for(auto j = 0u; j < tab.nrow; j++){
		auto time_str = tab.ele[j][0];
		auto t = number(time_str);
		if(t == UNSET){
			alert_source("The time '"+time_str+"' is not a number",so,0,j);
			return;
		}
		
		auto val_str = tab.ele[j][1];
		auto value = number(val_str);
		if(value == UNSET){
			alert_source("The value '"+val_str+"' is not a number",so,1,j);
			return;
		}
		
		auto col = 2;
		
		auto cf = so.comp_filt;
		for(auto cl = 0u; cl < ncla; cl++){
			const auto &claa = cla[cl];
			
			switch(cf.cla[cl].type){
			case ALL_FILT:
				for(auto c = 0u; c < claa.ncomp; c++){
					cf.cla[cl].comp[c] = true;
				}
				break;
			
			case FILE_FILT:
				col++;
				emsg("TO DO FILT");
				break;
			
			case COMP_FILT:
				break;
			}
		}	
	
		auto sd = 1.0;
		switch(so.obs_model.type){
		case PERCENT_OBSMOD: 
			sd = value*so.obs_model.percent/100; 
			break;
			
		case SD_OBSMOD: 
			sd = so.obs_model.sd;
			break;
			
		case FILE_OBSMOD: 
			{
				auto sd_str = tab.ele[j][col];
				auto sd = number(sd_str);
				if(sd == UNSET){
					alert_source("The standard deviation '"+sd_str+"' is not a number",so,2,j);
					return;	
				}
			}
			break;
		}
		if(sd < 1) sd = 1;
		
		PopData pd;
		pd.so = so.index;
		pd.t = t;
		pd.value = value;
		pd.sd = sd;
		pd.filt = global_convert(cf);
		pop_data.push_back(pd);
	}		
}


/// pop-trans-data command
void Species::population_trans_data(const DataSource &so)
{
	const auto &tab = so.table;
	
	for(auto j = 0u; j < tab.nrow-1; j++){
		auto time_str = tab.ele[j][0];
		auto t = number(time_str);
		if(t == UNSET){
			alert_source("The time '"+time_str+"' is not a number",so,0,j);
			return;
		}
		
		auto time_str_next = tab.ele[j+1][0];
		auto t_next = number(time_str_next);
		if(t_next == UNSET){
			alert_source("The time '"+time_str_next+"' is not a number",so,0,j+1);
			return;
		}
		
		auto val_str = tab.ele[j][1];
		auto value = number(val_str);
		if(value == UNSET){
			alert_source("The value '"+val_str+"' is not a number",so,1,j);
			return;
		}
		
		auto col = 2;
		
		auto tf = so.trans_filt;
		auto cf = so.comp_filt;
		for(auto cl = 0u; cl < ncla; cl++){
			const auto &claa = cla[cl];
				
			if(cl == so.cl){
				for(auto c = 0u; c < claa.ncomp; c++){
					cf.cla[cl].comp[c] = true;
				}
			}
			else{
				switch(cf.cla[cl].type){
				case ALL_FILT:
					for(auto c = 0u; c < claa.ncomp; c++){
						cf.cla[cl].comp[c] = true;
					}
					break;
				
				case FILE_FILT:
					col++;
					emsg("TO DO FILT");
					break;
				
				case COMP_FILT:
					break;
				}
			}	
		}
		
		auto sd = 1.0;
		switch(so.obs_model.type){
		case PERCENT_OBSMOD: 
			sd = value*so.obs_model.percent/100; 
			break;
			
		case SD_OBSMOD: 
			sd = so.obs_model.sd;
			break;
			
		case FILE_OBSMOD: 
			{
				auto sd_str = tab.ele[j][col];
				auto sd = number(sd_str);
				if(sd == UNSET){
					alert_source("The standard deviation '"+sd_str+"' is not a number",so,2,j);
					return;	
				}
			}
			break;
		}
		if(sd < 1) sd = 1;
		
		auto filt = trans_global_convert(so.cl,tf,cf);
		
		auto ti_min = get_ti(t);
		auto ti_max = get_ti(t_next);
		
		auto ref = pop_trans_data.size();
		for(auto ti = ti_min; ti < ti_max; ti++){
			for(auto tr = 0u; tr < tra_gl.size(); tr++){
				if(filt[tr] == true) pop_trans_ref[ti][tr].push_back(ref);
			}
		}
		
		PopTransData ptd;
		ptd.so = so.index;
		ptd.tmin = t;
		ptd.tmax = t_next;		
		ptd.value = value;
		ptd.sd = sd;
		ptd.filt = filt;
		pop_trans_data.push_back(ptd);
	}		
}


/// Finds the compartment from the name
unsigned int Species::find_c(unsigned int cl, string name) const
{
	const auto &claa = cla[cl];
	auto c = 0u; while(c < claa.ncomp && claa.comp[c].name != name) c++;
	if(c == claa.ncomp) return UNSET;
	return c;
}


/// Finds an individual with a given name, overwise creates a new individual
unsigned int Species::find_individual(string name)
{
	auto i = 0u; while(i < individual.size() && individual[i].name != name) i++;
	if(i == individual.size()){
		IndData ind;
		ind.name = name;
		individual.push_back(ind);
	}
	return i;
}


/// Used to generate some hypotheical data
void Species::generate_data() const
{
	
	/*
	ofstream fout("ind_data.csv");
	
	fout << "id,cinit,ds,sex,fe" << endl;
	for(auto i = 0u; i < 100; i++){
		fout << "Ind. " << i << ",";
		if(i < 98) fout << "S"; else fout << "I";
		fout << ",";
		if(ran() < 0.5) fout << "M"; else fout << "F";
		fout << ",";
		fout << normal_sample(0,1);
		fout << endl;
	}		
	cout << "DATA GENERATED" << endl;
	*/
}


/// Alerts a problem with a data source
void Species::alert_source(string st, const DataSource &so, unsigned int c, unsigned int r)
{
	if(c != UNSET) st += "(col '"+so.table.heading[c]+"', row "+tstr(r+1)+")";
	WarnData wa; wa.te = st; wa.line_num = so.line_num;
	warn.push_back(wa);
}


/// Orders X_vector such that order agrees with individual
void Species::X_vector_order()
{
	for(auto &fe : fix_effect){
		const auto &xvec_old = fe.X_vector;
		
		Xvector xvec_new;
		for(auto i = 0u; i < individual.size(); i++){
			auto name = individual[i].name;
			xvec_new.ind_list.push_back(name);
			
			auto j = find_in(xvec_old.ind_list,name,i);
			if(j == UNSET) xvec_new.value.push_back(UNSET);
			else xvec_new.value.push_back(xvec_old.value[j]);			
		}
		
		fe.X_vector = xvec_new;		
	}
}


/// Converts from a filter in 
vector <bool> Species::global_convert(const Filter &filt) const
{
	vector <bool> gfilt;
	
	for(auto c = 0u; c < comp_gl.size(); c++){
		const auto &cgl = comp_gl[c];
		
		auto cl = 0u; while(cl < ncla && filt.cla[cl].comp[cgl.cla_comp[cl]] == true) cl++;
		
		if(cl == ncla) gfilt.push_back(true);
		else gfilt.push_back(false);
	}
	
	return gfilt;
}


/// Converts from transition filter in classification to global filter
vector <bool> Species::trans_global_convert(unsigned int cl, const vector <bool> &trans_filt, const Filter &comp_filt)
{
	vector <bool> gfilt;
	
	auto filtg = global_convert(comp_filt); 
	for(auto tr = 0u; tr < tra_gl.size(); tr++){
		const auto &trg = tra_gl[tr];
		auto val = false;
		if(trg.cl == cl && trans_filt[trg.tr] == true){
			auto c = trg.i; if(c == SOURCE) c = trg.f;
			if(filtg[c] == true) val = true;
		}
		
		gfilt.push_back(val);
	}	
	
	return gfilt;
}


// Used to order events
bool EventData_ord (EventData ev1, EventData ev2)                      
{ return (ev1.t < ev2.t); };  

// Used to order observations
bool ObsData_ord (ObsData ob1, ObsData ob2)                      
{ return (ob1.t < ob2.t); };  

bool PopData_ord (PopData pd1, PopData pd2)                      
{ return (pd1.t < pd2.t); }; 

bool PopTransData_ord (PopTransData pd1, PopTransData pd2)                      
{ return (pd1.tmin < pd2.tmin); }; 

/// Orders individual events by time
void Species::order_data_events()
{
	sort(pop_data.begin(),pop_data.end(),PopData_ord);    
	sort(pop_trans_data.begin(),pop_trans_data.end(),PopTransData_ord);    
	
	// Adds populations onto the timelines of each individual
	for(auto j = 0u; j < pop_data.size(); j++){
		const auto &pd = pop_data[j];
		for(auto i = 0u; i < individual.size(); i++){
			auto &ind = individual[i];
			ObsData ob; 
			ob.type = OBS_POP;
			ob.ref = j;
			ob.t = pd.t;
			ind.obs.push_back(ob);
		}
	}
	
	
	// Sorts events by time
	for(auto i = 0u; i < individual.size(); i++){
		auto &ind = individual[i];
		sort(ind.ev.begin(),ind.ev.end(),EventData_ord);    
		sort(ind.obs.begin(),ind.obs.end(),ObsData_ord);    
	}
	
	if(false){
		for(auto i = 0u; i < individual.size(); i++){
			const auto &ind = individual[i];
			//print_event_data(ind.name,ind.ev);
			print_obs_data(ind.name,ind.obs);
		}
		
		
		for(const auto &pd : pop_data){
			cout << pd.t << " " << pd.value << " " << pd.sd << " Population data" << endl;
		}
		
		for(const auto &ptd : pop_trans_data){
			cout << ptd.tmin << " " << ptd.tmax << " " << ptd.value << " " << ptd.sd << " Population trans data" << endl;
		}
	}
}
