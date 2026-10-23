/// Routines related to validation of BICI (tornado plots and parameter scans)

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath> 
#include <algorithm> 
 
using namespace std;

#include "validation.hh"
#include "simulate.hh"
#include "data_sim.hh"
#include "utils.hh"

#ifndef MAC
#include <filesystem>
#endif

Validation::Validation()
{
}


/// Loads up a specification file
SpecInfo Validation::load_spec_info(SpecInfoType type, string spec_file)
{
	ifstream fin(spec_file);
	if(!fin) alert_input("File '"+spec_file+"' could not be loaded");
	
	vector <string> lines;
	
	string line;
	do{
		getline(fin,line);
		remove_cr(line);
		line = trim(line);
		if(line != "") lines.push_back(line);
		if(fin.eof()) break;
	}while(true);
		
	SpecInfo si;
	si.type = type;
	si.n = 20;
	si.log_on = false;
	si.data_param = false;
	si.para = true;
				
	double min = UNSET, max = UNSET;
	string wa = "In the spcification ";
		
	for(auto li : lines){
		if(begin_str(li,"{") && end_str(li,"}")){
			auto info = li.substr(1,li.length()-2);
			
			wa = "In specification '"+info+"' ";
			
			auto spl = split(info,' ');
			
			
			for(auto i = 1u; i < spl.size(); i++){
				auto va = trim(spl[i]);
				if(va != ""){
					auto spl2 = split(va,'=');
					
					if(spl2.size() != 2){
						if(spl2.size() == 1) run_error(wa+"syntax error in '"+va+"'. This should contain an '=' sign.");
						else run_error(wa+"syntax error in '"+va+"'. This should not contain multiple '=' signs.");
					}
					
					auto ta = spl2[0];
					auto value = trim(remove_quote(spl2[1]));

					auto done = false;
					
					if(ta == "param"){
						done = true;
						if(type != SCAN_PLOT) run_error(wa+"the 'param' option should only be specified for scan plots");
						
						if(begin_str(value,"$")){
							value = value.substr(1);
							si.data_param = true;
						}
						si.param_name = value;
					}
					
					if(ta == "range"){
						done = true;
						if(type != SCAN_PLOT) run_error(wa+"the 'range' option should only be specified for scan plots");
						
						auto spl3 = split(value,',');
						
						auto fl = false;
						if(spl3.size() != 2) fl = true;
						else{
							min = number(spl3[0]);
							max = number(spl3[1]);
							if(min == UNSET || max == UNSET) fl = true;
						}
						
						if(fl){
							run_error(wa+"the value for range '"+value+"' is not valid.");
						}
					}
					
					if(ta == "log"){
						done = true;
						if(value == "true") si.log_on = true;
						else{
							if(value == "false") si.log_on = false;
							else{
								run_error(wa+"the value for log '"+value+"' is not valid.");
							}
						}
					}
				
					if(ta == "n"){
						done = true;
						si.n = number(value);
						if(si.n == UNSET) run_error(wa+"the value for n '"+value+"' is not valid.");
						if(si.n <= 0) run_error(wa+"the value for n '"+value+"' must be positive.");
						if((unsigned int) si.n != si.n) run_error(wa+"the value for n '"+value+"' must be an integer.");
					}
					
					if(ta == "para"){
						done = true;
						if(value == "true") si.para = true;
						else{
							if(value == "false") si.para = false;
							else{
								run_error(wa+"the value for para '"+value+"' is not valid.");
							}
						}
					}	
							
					if(ta == "stop-extinction"){
						done = true;
						si.extinct_info = load_extinction(value);
					}
					
					if(!done){
						run_error(wa+"the value '"+spl2[0]+"' is not understood.");
					}
				}
			}
		}
		else{
			si.data_sim.push_back(li);
		}
	}

	if(type == SCAN_PLOT){
		if(si.param_name == "") run_error(wa+"a value for 'param' must be set.");
		
		if(min == UNSET || max == UNSET){
			run_error(wa+"a value for 'range' must be set.");
		}
	
		if(min >= max) run_error(wa+"the minimum must be smaller than the maximum.");
	
		if(si.log_on){
			if(min <= 0) run_error(wa+"the minimum must be positive when log is used.");
		}

		for(auto i = 0u; i <= si.n; i++){
			if(si.log_on){
				si.value_sim.push_back(exp(log(min)+(log(max)-log(min))*double(i)/si.n));	
			}
			else{
				si.value_sim.push_back(min+(max-min)*double(i)/si.n);
			}
		}
	}
	
	return si;
}


/// Deletes any exisiting files or observation parameters
void Validation::delete_files(const Model &model, Output &output, string tor_dir, string run_dir, string op_dir, string spec_file)
{
	// Removes any observation model parameters
	for(const auto &par : model.param){
		if(par.mod_param == OBS_MOD_PARAM){
			output.delete_command(par.name,par.line_num,false);
		}
	}

	// Removes any existing data 
	for(auto p = 0u; p < model.nspecies; p++){
		for(const auto &so : model.species[p].source){
			output.delete_command(so.name,so.line_num,false);
		}
	}
	
	
	output.ensure_directory(tor_dir);
	output.ensure_directory(run_dir);
	output.ensure_directory(op_dir);
	
	// Copies spec_file into output directory
	ifstream fin(spec_file);
	if(!fin) alert_input("File '"+spec_file+"' could not be loaded");
	
	ofstream fout(tor_dir+"/spec_file.txt");

	string line;
	do{
		getline(fin,line);
		remove_cr(line);
		line = trim(line);
		fout << line << endl;
		if(fin.eof()) break;
	}while(true);
}

		
/// Sets up a tornado run
void Validation::tornado_setup(string spec_file, ExtFactor ext_factor, string file) 
{
	auto nc = num_core();
	if(nc != 1) run_error("Setup requires only one core");
	
	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	
	auto name = file.substr(0,file.length()-5);
	
	auto tor_dir = name+"_tornado/";
	auto run_dir = tor_dir+"run/";
	auto op_dir = tor_dir+"op/";
	
	auto ncore_inf = 1u;
	
	auto si = load_spec_info(TORNADO_PLOT,spec_file);
		
	auto valid_run = get_valid_run(si.n,run_dir,si.data_sim);

	{                                        // Removes data and copies into 
		Model model(INF,ext_factor,true); 
		Mpi mpi(UNSET,model);   
		Output output(model,mpi);  
	
		Input input(model,file,UNSET,mpi);     // Imports information from input into model

		if(input.datadir != "") alert_input("Tornado can currently only be run on BICI-script with embedded files (i.e. no separate data directory).");
	
		output.init(input);

		const auto &det = model.details;
		
		ncore_inf = det.nchain/det.num_per_core;
			
		delete_files(model,output,tor_dir,run_dir,op_dir,spec_file);
		
		for(auto r = 0u; r < valid_run.size(); r++){
			output.copy(valid_run[r].file);
		}
	}
	
	run_sim_data(name,op_dir,si.para,valid_run,ncore_inf,ext_factor,file,"tornado-res",spec_file,si.extinct_info);
}


/// Gets a list of files for replicates
vector <ValidRun> Validation::get_valid_run(unsigned int n, string run_dir, const vector <string> &data_sim_lines) const
{
	vector <ValidRun> valid_run;
	for(auto r = 0u; r < n; r++){
		ValidRun vr; 
		vr.file = run_dir+"Run_"+tstr(r+1)+".bici";
		vr.data_sim_lines = data_sim_lines;
		valid_run.push_back(vr);
	}
	
	return valid_run;
}

	
/// Runs simulations and generates data
void Validation::run_sim_data(string name, string op_dir, bool para_run, const vector <ValidRun> &valid_run, unsigned int ncore_inf, ExtFactor ext_factor, string file, string com, string spec_file, const ExtinctInfo &extinct_info) const
{
	percentage_start(RUN_PER);
	for(auto r = 0u; r < valid_run.size(); r++){  
		percentage(r,valid_run.size());
		
		const auto &vr = valid_run[r]; 
		
		auto fi = vr.file;
		
		{   // Simulates state
			Model model(SIM,ext_factor,true); 
			Mpi mpi(UNSET,model);   
			Output output(model,mpi,true);    
			Input input(model,fi,r,mpi,true);  
			output.init(input);
	
			Simulate simu(model,output,mpi,true);
			auto loop = 0u, loop_max = 100u;
			do{
				simu.run();
				if(!simu.extinction_criterion(extinct_info)) break;
				output.clear();
				loop++;
			}while(loop < loop_max);
			if(loop == loop_max) run_error("Could not satisfy 'stop-extinction' criterion within "+tstr(loop_max)+" simulations from the model.");
				
			output.end(fi,UNSET);
		}
		
		{ // Generates data
			Model model(DATA_SIM,ext_factor,true); 
			Mpi mpi(UNSET,model);   
			Output output(model,mpi,true);    
			Input input(model,fi,r,mpi,true);  
			output.init(input);
			
			DataSim data_sim(model,output,true);
			bool replace_data = false;
			for(auto &ds : vr.data_sim_lines){
				data_sim.run(ds,1,replace_data);
			}
			output.end(fi,UNSET);
		}
	}
	percentage_end();
	
	auto spl = split(name,'/');
	
	auto scr = spl[spl.size()-1]+"_script";             // Output the script for running inference
	scr += ".sh";
	
	string code = "./bici-para";
#ifdef WINDOWS
	code = "./bici-core.exe";
#endif
	
	ofstream fout(scr);		
	check_open(fout,scr);
	
	for(auto r = 0u; r < valid_run.size(); r++){
		auto fi = valid_run[r].file; 
		
		fout << "echo ''" << endl;
		fout << "echo ''" << endl;
		fout << "echo 'Run " << fi << "'" << endl;
		if(para_run) fout << "nohup ";
		if(ncore_inf > 1) fout << "mpirun -n " << ncore_inf << " "; 

		fout << code << " " << fi << " -seed=" << r << " inf";
		if(para_run) fout << " > " << op_dir << "inf_" << r << ".txt 2>&1&";
		fout << endl << endl;
	}
	
	fout << "echo ''" << endl;
	
	if(para_run) fout << "echo \"WHEN INFERENCE COMPLETE: " << code << " " << file << " " << com << " " << spec_file << "\"" << endl; 
	else fout << "echo \"COMBINE RESULTS: " << code << " " << file << " " << com << " " << spec_file << "\"" << endl; 
		
	cout << endl;

#ifdef WINDOWS
#else
	cout << "MAKE EXECUTABLE: chmod +x " << scr << endl;
#endif

	cout << "RUN SCRIPT: ./" << scr << endl; 
}

	
/// Used to order proposal speeds
bool Stat_ord (const Stat &st1, const Stat &st2)                      
{ return (st1.mean < st2.mean); };  


/// Makes sure each table has the same row (otherwise missing data could cause different sizes) 
vector <Table> Validation::consistent_row(const vector <Table> &list) const
{
	vector <string> row_sel;
	
	for(const auto &ta : list){
		if(ta.nrow > row_sel.size()){
			row_sel.clear();
			for(const auto &row : ta.ele) row_sel.push_back(row[0]);
		}
	}
	
	auto N = row_sel.size();
	
	vector <Table> list_new;
	
	for(const auto &ta : list){
		auto ta_new = ta;
		ta_new.ele.clear();
		ta_new.ele.resize(N);
		
		for(auto r = 0u; r < ta.nrow; r++){
			auto j = find_in(row_sel,ta.ele[r][0]);
			if(j == UNSET) emsg("Cannot find row"); 
			ta_new.ele[j] = ta.ele[r];
		}
		
		for(auto r = 0u; r < N; r++){
			if(ta_new.ele[r].size() == 0){
				ta_new.ele[r].push_back(row_sel[r]);
				for(auto c = 1u; c < ta.ncol; c++) ta_new.ele[r].push_back("NA");
			}
		}
		
		list_new.push_back(ta_new);
	}
	
	return list_new;
}


/// Accounts for unset value
string Validation::us(double val) const
{
	if(val == UNSET) return "NA";
	return tstr(val);
}	
	

/// Combines results from a tornado run
void Validation::tornado_result(string spec_file, ExtFactor ext_factor, bool no_question, string file)
{
	auto si = load_spec_info(TORNADO_PLOT,spec_file);
	
	Model model_base(SIM,ext_factor,no_question);  
	Mpi mpi_base(UNSET,model_base);
	Output op(model_base,mpi_base);      
	Input input(model_base,file,UNSET,mpi_base,true);
	
	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	auto name = file.substr(0,file.length()-5);
	
	auto tor_dir = name+"_tornado/";
	auto run_dir = tor_dir+"run/";
	auto op_dir = tor_dir+"op/";
	
	if(!op.directory_exist(tor_dir)){
		run_error("The directory '"+tor_dir+"' does not exist. Has analysis been performed?");
	}
	
	if(!op.directory_exist(run_dir)){
		run_error("The directory '"+run_dir+"' does not exist. Has analysis been performed?");
	}
	
	auto nrun = si.n;
	
	auto valid_run = get_valid_run(nrun,run_dir,vector<string>());
		
	vector <Table> param_stats;
	vector <Table> pred_acc_stats;
		
	percentage_start(RUN_PER);
	for(auto r = 0u; r < valid_run.size(); r++){
		percentage(r,valid_run.size());
	
		auto fi = valid_run[r].file;
		
		Model model(INF,ext_factor,no_question);  
		Mpi mpi(UNSET,model);
		Output output(model,mpi,true);           
		Input input(model,fi,UNSET,mpi,true);    
		
		{
			const auto &tab = model.inf_param_stats;
			if(tab.heading.size() == 0) run_error("Parameter statistics could not be found. Have you run the script?");
		
			param_stats.push_back(tab);
		}
		
		{
			const auto &tab = model.inf_pred_acc;
			if(tab.heading.size() > 0){
				pred_acc_stats.push_back(tab);
			}
		}
	}
	
	param_stats = consistent_row(param_stats);
	
	percentage_end();
	
	{
		auto R = param_stats.size(); 
		
		const auto &tab0 = param_stats[0];
		auto N = tab0.nrow;
		
		vector <bool> value_unset(N,false);

		vector <double> value;
		{ // Works out the true value for the parameter
			for(auto i = 0u; i < N; i++){
				auto na = tab0.ele[i][0];
				auto val = get_param_value(na,model_base);
				value.push_back(val);
			}
		}
	
		vector < vector <Stat> > stat;
		stat.resize(N);
		
		{	
			auto file = tor_dir+"ESS.csv";
			ofstream fout(file);
			check_open(fout,file);
			
			fout << "Run";

			for(auto i = 0u; i < N; i++) fout << ",\"" << add_escape_char(tab0.ele[i][0]) << " ESS\"";
			fout << endl;
			
			for(auto r = 0u; r < R; r++){
				const auto &tab = param_stats[r];
				fout << (r+1);
				for(auto i = 0u; i < N; i++){
					Stat st; 
					st.mean = number(tab.ele[i][1]);
					st.sd = number(tab.ele[i][2]);
					st.CImin = number(tab.ele[i][3]);
					st.CImax = number(tab.ele[i][4]);
					stat[i].push_back(st);
					
					fout << "," << tab.ele[i][5];
				}				
				fout << endl;
			}
		}
		
		{	
			auto file = tor_dir+"GR.csv";
			ofstream fout(file);
			check_open(fout,file);
			
			fout << "Run";
			for(auto i = 0u; i < N; i++) fout << ",\"" << add_escape_char(tab0.ele[i][0]) << " GR\"";
			fout << endl;
			
			for(auto r = 0u; r < R; r++){
				const auto &tab = param_stats[r];
				
				fout << (r+1);
				for(auto i = 0u; i < N; i++) fout << "," << tab.ele[i][6];
				fout << endl;
			}
		}
		
		auto tname = "result.csv";
		{	
			auto file = tor_dir+tname;
			ofstream fout(file);
			check_open(fout,file);
			
			fout << "Sorted Run";
			for(auto i = 0u; i < N; i++){
				auto na = add_escape_char(tab0.ele[i][0]);
				fout << ",\"" << na << " value\"" << ",\"" << na << " (post. mean)\"" << ",\"" << na << " (CI min)\"" << ",\"" << na << " (CI max)\"";
			}
			fout << endl;
			
			for(auto i = 0u; i < N; i++){
				sort(stat[i].begin(),stat[i].end(),Stat_ord);
			}
				
			for(auto r = 0u; r < R; r++){
				fout << (r+1);
				for(auto i = 0u; i < N; i++){
					fout << ",";
					if(value[i] == UNSET) value_unset[i] = true;
					fout << us(value[i]) << "," << us(stat[i][r].mean) << "," << us(stat[i][r].CImin) << "," << us(stat[i][r].CImax);
				}
				fout << endl;
			}
		}
		
		auto un_name = "uncertainty.csv";
		{	
			auto file = tor_dir+un_name;
			ofstream fout(file);
			check_open(fout,file);
			
			for(auto i = 0u; i < N; i++){
				auto na = add_escape_char(tab0.ele[i][0]);
				if(i != 0) fout << ",";
				fout << "\"" << na << "\"";
			}
			fout << endl;
		
			for(auto i = 0u; i < N; i++){
				if(i != 0) fout << ",";
				
				auto sd_av = 0.0, mean_av = 0.0;
				for(auto r = 0u; r < R; r++){
					mean_av += stat[i][r].mean;
					sd_av += stat[i][r].sd;
				}
				mean_av /= R;
				sd_av /= R;
			
				fout << sd_av/mean_av;
			}		
			fout << endl;
		}
		
		auto gp = tor_dir+"gp/";
		op.ensure_directory(gp);

		{ // Creates lines for true value
			auto file = gp+"true_value.csv";
			ofstream fout(file);
			check_open(fout,file);
			
			fout << "y";
			for(auto i = 0u; i < N; i++){
				auto na = add_escape_char(tab0.ele[i][0]);
				fout << ",\"" << na << " value\"";
			}
			fout << endl;
			fout << 0;
			for(auto i = 0u; i < N; i++) fout << "," << us(value[i]);
			fout << endl;
			fout << nrun+1;
			for(auto i = 0u; i < N; i++) fout << "," << us(value[i]);
			fout << endl;
		}
		
		{  // Creates gp scripts
			auto file = gp+"gp_single.txt";
			ofstream fout(file);
			check_open(fout,file);
		
			fout << gp_header();
				
			fout << "set yrange [0:" << nrun+1 << "]" << endl;
			fout << "unset ytics" << endl;

			for(auto i = 0u; i < N; i++){
				auto name = add_escape_char(tab0.ele[i][0]);
				
				name = convert_label(name);
				
				fout << "set xlabel '"+name+"' font ',30'" << endl;
				fout << "unset ylabel" << endl;

				fout << "plot '../" << tname << "' using " << 3+4*i << ":1 ps 1 lw 4 lc '#000000' notitle,\\" << endl;
				fout << "'../" << tname << "' using " << 3+4*i << ":1:" << 3+4*i+1 << ":" << 3+4*i+2 << " with xerrorbars ls 1 lw 4 notitle,\\" << endl;
				fout << "'true_value.csv' using " << 2+i << ":1 with lines ls 2 notitle" << endl;
			}
		}
		
		{  // Creates gp scripts
			auto file = gp+"gp.txt";
			ofstream fout(file);
			check_open(fout,file);
		
			fout << gp_header_multi();
				
			fout << "set yrange [0:" << nrun+1 << "]" << endl;
			fout << "unset ytics" << endl;

			for(auto i = 0u; i < N; i++){
				auto name = add_escape_char(tab0.ele[i][0]);
				
				name = convert_label(name);
				
				fout << "set xlabel '"+name+"' font ',15'" << endl;
				fout << "unset ylabel" << endl;

				fout << "plot '../" << tname << "' using " << 3+4*i << ":1 ps 1 lw 2 lc '#000000' notitle,\\" << endl;
				fout << "'../" << tname << "' using " << 3+4*i << ":1:" << 3+4*i+1 << ":" << 3+4*i+2 << " with xerrorbars ls 1 lw 2 notitle,\\" << endl;
				fout << "'true_value.csv' using " << 2+i << ":1 with lines ls 2 lw 2 notitle" << endl;
			}
		}
	}
	
	{
		auto R = pred_acc_stats.size(); 
		
		auto tab0 = pred_acc_stats[0];
		for(auto r = 0u; r < tab0.nrow; r++){
			vector <double> vec;
			for(auto j = 0u; j < R; j++){
				auto val = number(pred_acc_stats[j].ele[r][1]);
				if(val == UNSET) emsg("Prediction accuracy not a number");
				vec.push_back(val);
			}
			
			auto stat = get_statistic(vec);
			
			stringstream ss;
			ss << stat.mean << " (" << stat.CImin << " - " << stat.CImax << ")";
			tab0.ele[r][1] = ss.str();
		}
		
		auto tname = "pred_acc.csv";
		
		auto file = tor_dir+tname;
		ofstream fout(file);
		check_open(fout,file);
		
		fout << output_table(tab0);
	}

	cout << endl << "RESULTS: In directory '" << tor_dir << "'." << endl;
	//cout << "PLOT: Gnuplot can be used to plot results. Go to the 'gp' sub-directory and type 'gnuplot gp.txt' to generate graphs and 'ps2pdf gp.ps gp.pdf' to convert to pdf." << endl;
	
	cout << "PLOT: Gnuplot can be used to plot results. Go to the 'gp' sub-directory and type 'gnuplot gp_single.txt' to generate graphs and 'ps2pdf gp_single.ps gp.pdf' to convert to pdf." << endl;
}


/// Sets up a parameter scan
void Validation::scan_setup(string spec_file, ExtFactor ext_factor, string file) 
{
	auto nc = num_core();
	if(nc != 1) run_error("Setup requires only one core");
	
	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	
	auto name = file.substr(0,file.length()-5);
	
	auto si = load_spec_info(SCAN_PLOT,spec_file);
		
	auto na = add_escape_char(si.param_name);
	na = replace(na,"\\","");

	auto tor_dir = name+"_scan_"+na+"/";
	auto run_dir = tor_dir+"run/";
	auto op_dir = tor_dir+"op/";
		
	auto ncore_inf = 1u;
	
	vector <ValidRun> valid_run;
	{                                        // Removes data and copies into 
		Model model(INF,ext_factor,true); 
		Mpi mpi(UNSET,model);   
		Output output(model,mpi);    
		Input input(model,file,UNSET,mpi);     // Imports information from input into model
		
		if(input.datadir != "") alert_input("Scanning can currently only be run on BICI-script with embedded files (i.e. no separate data directory).");
		
		output.init(input);

		if(!si.data_param){                    // Checks parameter exists and makes sure values within prior
			auto pr = get_param_ref(si.param_name,model);
			if(pr.th == UNSET || pr.index == UNSET){
				run_error("Parameter '"+si.param_name+"' could not be found in the model.");
			}
			
			auto pb = model.param_bounds(pr);
			
			if(pb.min != UNSET && pb.max != UNSET){
				for(auto va : si.value_sim){
					if(va <= pb.min || va >= pb.max){
						run_error("Parameter values for '"+si.param_name+"' are outside of the range ("+tstr(pb.min)+" - "+tstr(pb.max)+") specified by the prior.");
					}
				}
			}
		}
		
		valid_run = get_valid_run(si.value_sim.size(),run_dir,si.data_sim);
		
		const auto &det = model.details;
	
		ncore_inf = det.nchain/det.num_per_core;
	
		delete_files(model,output,tor_dir,run_dir,op_dir,spec_file);
	
		for(auto r = 0u; r < valid_run.size(); r++){
			if(si.data_param){ // If parameter in data specification then alters that
				auto sub = false;
				auto rep = tstr(si.value_sim[r]);
				auto na = "$"+si.param_name;
				auto len = na.length();
				for(auto &st : valid_run[r].data_sim_lines){			
					auto i = 0u;
					while(i+len < st.length()){
						if(st.substr(i,len) == na){
							sub = true;
							st = st.substr(0,i)+rep+st.substr(i+len);
						}
						i++;
					}
				}
				
				if(!sub) run_error("Could not find parameter '"+na+"' in the data specification.");
			}				
			else{              // If parameter in BICI-script then alter that
				output.change_sim_value(si.param_name,si.value_sim[r]);
			}
			output.copy(valid_run[r].file);
		}
	}
	
	run_sim_data(name+"_"+na,op_dir,si.para,valid_run,ncore_inf,ext_factor,file,"scan-res",spec_file,si.extinct_info);
}


/// Makes sure not all values are zero (otherwise gnuplot interprets as undefined)
bool Validation::is_log(const SpecInfo &si, vector <Stat> &stat) const
{
	if(si.log_on == false) return false;
	
	auto fl = false;
	for(auto r = 0u; r < stat.size(); r++){
		if(stat[r].mean <= 0) fl = true;
		if(stat[r].CImin <= 0) fl = true;
		if(stat[r].CImax <= 0) fl = true;
	}
	
	if(fl == true) return false;
	return true;
}


/// Determines label on x axis
string Validation::get_si_name(const SpecInfo &si) const
{
	if(si.data_param) return convert_label(si.param_name);
	return convert_label("True "+si.param_name);
}

		
/// Combines results from scan runs
void Validation::scan_result(string spec_file, ExtFactor ext_factor, bool no_question, string file)
{
	auto si = load_spec_info(SCAN_PLOT,spec_file);
	
	Model model_base(SIM,ext_factor,no_question);  
	Mpi mpi_base(UNSET,model_base);
	Output op(model_base,mpi_base);      
	Input input(model_base,file,UNSET,mpi_base,true);

	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	auto name = file.substr(0,file.length()-5);
	
	auto pname = si.param_name;
	
	auto na = add_escape_char(pname);
	na = replace(na,"\\","");
	auto tor_dir = name+"_scan_"+na+"/";
	auto run_dir = tor_dir+"run/";
	
	if(!op.directory_exist(tor_dir)){
		run_error("The directory '"+tor_dir+"' does not exist. Has analysis been performed?");
	}
	
	if(!op.directory_exist(run_dir)){
		run_error("The directory '"+run_dir+"' does not exist. Has analysis been performed?");
	}

	auto prop = get_param_prop(si.param_name);
	
	auto valid_run = get_valid_run(si.value_sim.size(),run_dir,vector<string>());
		
	vector <Table> param_stats;
	
	percentage_start(RUN_PER);
	for(auto r = 0u; r < valid_run.size(); r++){
		percentage(r,valid_run.size());
	
		auto fi = valid_run[r].file;
		
		Model model(INF,ext_factor,no_question);  
		Mpi mpi(UNSET,model);
		Output output(model,mpi,true);           
		Input input(model,fi,UNSET,mpi,true);    
		
		const auto &tab = model.inf_param_stats;
		if(tab.heading.size() == 0) run_error("Parameter statistics could not be found. Have you run the script?");
		
		param_stats.push_back(tab);
	}
	
	param_stats = consistent_row(param_stats);
		
	percentage_end();
	
	auto R = param_stats.size(); 
	
	const auto &tab0 = param_stats[0];
	auto N = tab0.nrow;
	
	vector <double> value;
	{ // Works out the true value for the parameter
		for(auto i = 0u; i < N; i++){
			auto na = tab0.ele[i][0];
			auto val = get_param_value(na,model_base);
			value.push_back(val);
		}
	}
	
	vector < vector <Stat> > stat;
	stat.resize(N);
	
	{	
		auto file = tor_dir+"ESS.csv";
		ofstream fout(file);
		check_open(fout,file);
		
		fout << "Run,\"" << pname << " value\"";
		for(auto i = 0u; i < N; i++) fout << ",\"" << add_escape_char(tab0.ele[i][0]) << " ESS\"";
		fout << endl;
		
		for(auto r = 0u; r < R; r++){
			const auto &tab = param_stats[r];
			
			fout << (r+1) << "," << si.value_sim[r];
			for(auto i = 0u; i < N; i++){
				Stat st; 
				st.mean = number(tab.ele[i][1]);
				st.sd = number(tab.ele[i][2]);
				st.CImin = number(tab.ele[i][3]);
				st.CImax = number(tab.ele[i][4]);
				stat[i].push_back(st);
				
				fout << "," << tab.ele[i][5];
			}				
			fout << endl;
		}
	}
	
	{	
		auto file = tor_dir+"GR.csv";
		ofstream fout(file);
		check_open(fout,file);
		
		fout << "Run,\"" << pname << " value\"";
		for(auto i = 0u; i < N; i++) fout << ",\"" << add_escape_char(tab0.ele[i][0]) << " GR\"";
		fout << endl;
		
		for(auto r = 0u; r < R; r++){
			const auto &tab = param_stats[r];
			
			fout << (r+1) << "," << si.value_sim[r];
			for(auto i = 0u; i < N; i++) fout << "," << tab.ele[i][6];
			fout << endl;
		}
	}

	vector <bool> value_unset(N,false);		

	auto tname = "result.csv";
	{	
		auto file = tor_dir+tname;
		ofstream fout(file);
		check_open(fout,file);
		
		fout << "Run," << si.param_name;
		for(auto i = 0u; i < N; i++){
			auto na = add_escape_char(tab0.ele[i][0]);
			auto prop2 = get_param_prop(na);
		
			fout << ",\"" << na << " value\"" << ",\"" << na << " (post. mean)\"" << ",\"" << na << " (CI min)\"" << ",\"" << na << " (CI max)\"";
		}
		fout << endl;

		for(auto r = 0u; r < R; r++){
			fout << (r+1);
			fout << "," << si.value_sim[r];
			for(auto i = 0u; i < N; i++){
				fout << ",";
				
				if(value[i] == UNSET)	value_unset[i] = true;
				fout << us(value[i]); 
			
				fout << "," << us(stat[i][r].mean) << "," << us(stat[i][r].CImin) << "," << us(stat[i][r].CImax);
			}
			fout << endl;
		}
	}

	auto gp = tor_dir+"gp/";
	op.ensure_directory(gp);
	
	{  // Creates gp scripts
		auto file = gp+"gp_single.txt";
		ofstream fout(file);
		check_open(fout,file);
	
		fout << gp_header();

		if(si.log_on) fout << "set logscale x" << endl; 
				
		fout << "set size square" << endl;
		
		auto namex = get_si_name(si);
			
		for(auto i = 0u; i < N; i++){
			auto na = tab0.ele[i][0];
			auto namey = add_escape_char(na);
			namey = convert_label(namey);
			
			if(!is_tt_stat(na) && is_log(si,stat[i])) fout << "set logscale y" << endl; 
			else fout << "unset logscale y" << endl; 
		
			fout << "set xlabel '"+namex+"' font ',30'" << endl;
			fout << "set ylabel '"+namey+"' font ',30'" << endl;
			
			fout << "plot '../" << tname << "' using " << 2 << ":" << 4+4*i << " ps 1 lw 4 lc '#000000' notitle,\\" << endl;
			fout << "'../" << tname << "' using " << 2 << ":" << 4+4*i << ":" << 4+4*i+1 << ":" << 4+4*i+2 << " with yerrorbars ls 1 lw 4 notitle";
			if(!value_unset[i]){			
				fout << ",\\" << endl;
				fout << "'../" << tname << "' using " << 2 << ":" << 3+4*i << " with lines ls 2 notitle";
			}
			fout << endl;
		}
	}

	{  // Creates gp scripts
		auto file = gp+"gp.txt";
		ofstream fout(file);
		check_open(fout,file);
	
		fout << gp_header_multi();
		if(si.log_on){
			fout << "set logscale x" << endl; 
			fout << "set logscale y" << endl; 
		}
		fout << "set size square" << endl;
		
		auto namex = get_si_name(si);
			
		for(auto i = 0u; i < N; i++){
			auto namey = add_escape_char(tab0.ele[i][0]);
			namey = convert_label(namey);
			
			fout << "set xlabel '"+namex+"' font ',15'" << endl;
			fout << "set ylabel '"+namey+"' font ',15'" << endl;
			
			fout << "plot '../" << tname << "' using " << 2 << ":" << 4+4*i << " ps 1 lw 2 lc '#000000' notitle,\\" << endl;
			fout << "'../" << tname << "' using " << 2 << ":" << 4+4*i << ":" << 4+4*i+1 << ":" << 4+4*i+2 << " with yerrorbars ls 1 lw 2 notitle,\\" << endl;
			fout << "'../" << tname << "' using " << 2 << ":" << 3+4*i << " with lines ls 2 lw 2 notitle" << endl;
		}
	}
	
	cout << endl << "RESULTS: In directory '" << tor_dir << "'." << endl;
	cout << "PLOT: Gnuplot can be used to plot results. Go to the 'gp' sub-directory and type 'gnuplot gp.txt' to generate graphs and 'ps2pdf gp.ps gp.pdf' to convert to pdf." << endl;
}


/// Returns a header for a gp file
string Validation::gp_header() const
{
	stringstream ss;
	ss << "set style line 1 lt 1 lc rgb '#000000' lw 4" << endl;
	ss << "set style line 2 lt 1 dt 3 lc rgb '#0000ff' lw 4" << endl;
	ss << "set datafile separator ','" << endl;
	ss << "set terminal postscript enhanced color font ',24' dl 3" << endl;
	ss << "set output 'gp_single.ps'" << endl;
	ss << "set bar 2" << endl;

	return ss.str();
}


/// Returns a header for a gp file
string Validation::gp_header_multi() const
{
	stringstream ss;
	ss << "set style line 1 lt 1 lc rgb '#000000' lw 2" << endl;
	ss << "set style line 2 lt 1 dt 3 lc rgb '#0000ff' lw 2" << endl;
	ss << "set datafile separator ','" << endl;
  ss << "set terminal postscript  portrait  enhanced color font ',12' dl 1.5" << endl;
	ss << "set output 'gp.ps'" << endl;
	ss << "set bar 1" << endl;
	ss << "set multiplot layout 3,2 rowsfirst margins 0.1,0.9,0.1,0.9 spacing 0.1" << endl;

	return ss.str();
}


/// Converts label into what can be used in gnuplot
string Validation::convert_label(string name) const 
{
	name = replace(name,"\\alpha","{/Symbol a}");  
	name = replace(name,"\\beta","{/Symbol b}");  
	name = replace(name,"\\gamma","{/Symbol g}"); 
	name = replace(name,"\\Gamma","{/Symbol G}"); 
	name = replace(name,"\\delta","{/Symbol d}"); 
	name = replace(name,"\\Delta","{/Symbol D}"); 
	name = replace(name,"\\epsilon","{/Symbol-Oblique e}"); 
	name = replace(name,"\\eta","{/Symbol-Oblique h}"); 
	name = replace(name,"\\Eta","{/Symbol-Oblique H}"); 
	name = replace(name,"\\theta","{/Symbol-Oblique q}");
	name = replace(name,"\\lambda","{/Symbol-Oblique l}");
	name = replace(name,"\\mu","{/Symbol m}"); 
	name = replace(name,"\\omega","{/Symbol-Oblique w}");
	name = replace(name,"\\pi","{/Symbol p}"); 
	name = replace(name,"\\sigma","{/Symbol s}"); 
	name = replace(name,"\\Omega","{/Symbol W}"); 
	name = replace(name,"\\phi","{/Symbol-Oblique f}");
	name = replace(name,"\\Phi","{/Symbol-Oblique F}");
	name = replace(name,"\\zeta","{/Symbol z}"); 
	
	name = replace(name,"\\Theta","{/Symbol-Oblique Q}");
	name = replace(name,"\\kappa","{/Symbol k}"); 
	name = replace(name,"\\Lambda","{/Symbol L}"); 
	name = replace(name,"\\iota","{/Symbol i}"); 
	name = replace(name,"\\nu","{/Symbol n}"); 
	name = replace(name,"\\xi","{/Symbol x}"); 
	name = replace(name,"\\Xi","{/Symbol X}"); 
	name = replace(name,"\\Pi","{/Symbol P}"); 
	name = replace(name,"\\rho","{/Symbol r}"); 
	name = replace(name,"\\tau","{/Symbol t}"); 
	name = replace(name,"\\omicron","{/Symbol o}"); 
	name = replace(name,"\\upsilon","{/Symbol u}"); 
	name = replace(name,"\\psi","{/Symbol y}"); 
	name = replace(name,"\\Psi","{/Symbol Y}"); 
	name = replace(name,"\\chi","{/Symbol c}"); 
	name = replace(name,"\\Chi","{/Symbol C}"); 
	
	auto spl = split(name,'^');
	if(spl.size() == 2){
		name = spl[0]+"^{"+spl[1]+"}";
	}
	return name;
}


/// Determines if a quantity is a statistic
bool Validation::is_tt_stat(string na) const 
{
	if(na == "N^origin" || na == "N^infected" || na == "N^unobs" || 
	   na == "N^mut-tree" || na == "N^mut-origin" || na == "t^root" || na == "α^IPA") return true;
	
	return false;
}


/// Gets tha value of a parameter
double Validation::get_param_value(string na, const Model &model) const 
{
	if(is_tt_stat(na)) return UNSET;
	
	auto pr = get_param_ref(na,model);
	
	const auto &par = model.param[pr.th];
	const auto &ef = par.element_ref[pr.index];
	if(ef.cons != true) run_error("Could not find value for parameter '"+na+"'");
	
	return model.constant.value[ef.index];
}

			
/// Gets a parameter reference from a parameter name
ParamRef Validation::get_param_ref(string na, const Model &model) const
{;
	auto prop = get_param_prop(na);
	const auto &param = model.param;
	
	auto pna = add_escape_char(prop.name);
	
	auto th = 0u; while(th < param.size() && add_escape_char(param[th].name) != pna) th++;
	if(th == param.size()) run_error("Could not find parameter '"+na+"'");
	
	auto index = 0u;
	
	const auto &par = param[th];
	if(par.N == 1){
		if(prop.dep.size() != 0) emsg("Problem loading parameter '"+na+"'");
	}
	else{
		if(prop.dep.size() != par.dep.size()) run_error("Dependency wrong for parameter '"+na+"'");
		
		for(auto d = 0u; d < par.dep.size(); d++){
			const auto &de = par.dep[d];
			auto j = find_in(de.list,prop.dep[d]);
			if(j == UNSET) emsg("Problem with dependency");
			index += de.mult*j;
		}
	}
	
	ParamRef pr; pr.th = th; pr.index = index;
	
	return pr;
}


/// Loads up information about stopping extinction
ExtinctInfo Validation::load_extinction(string te) const
{
	te =  trim(remove_quote(te));
	
	auto spl = split(te,',');
	
	if(spl.size() > 2) alert_input("The 'stop-extinction' definition '"+te+"' is not recognised");
	
	auto tra = spl[0];

	auto spl2 = split(tra,'>');
	if(spl2.size() != 2) alert_input("The 'stop-extinction' transition '"+tra+"' is not recognised");
	if(!end_str(spl2[0],"-")) alert_input("The 'stop-extinction' transition '"+tra+"' is not recognised");	
	
	ExtinctInfo ei;
	ei.on = true;
	ei.trans = tra;
	ei.comp_i = trim(spl2[0].substr(0,spl2[0].size()-1));
	ei.comp_f = trim(spl2[1]);

	if(spl.size() == 2){
		auto num_str = spl[1];
		auto num = number(num_str);
		if(num == UNSET) alert_input("The 'stop-extinction' transition number '"+num_str+"' is not recognised");
		if(num <= 0) alert_input("The 'stop-extinction' transition number '"+num_str+"' must be positive");
		if(num != (int) num) alert_input("The 'stop-extinction' transition number '"+num_str+"' must be an integer");
		
		ei.min_trans = num;
	}
	else{
		ei.min_trans = 4;
	}

	return ei;
}


