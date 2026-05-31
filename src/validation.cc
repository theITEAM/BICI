/// Routines related to validation of BICI (tornado plots etc)

// <<< Population-based model >>>
//./bici-para Execute/init.bici tornado 'ic' 'pop-data name="pop" filter="DS=file" dt=10 error="poisson"'
// ./bici-para Execute/init.bici scan:'\beta range=0.0005,0.002 n=20 log=false'  'ic' 'pop-data name="pop" filter="DS=file" dt=10 error="poisson"'


// <<< Individual-based model >>> 
//./bici-para Execute/init.bici tornado 'ic' 'pop-data name="pop" filter="DS=file" dt=10 error="poisson"'
// ./bici-para Execute/init.bici scan:'\beta'  'ic' 'pop-data name="pop" filter="DS=file" dt=10 error="poisson"'



#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath> 
#include <algorithm> 
#include <filesystem>
 
using namespace std;

#include "validation.hh"
#include "simulate.hh"
#include "data_sim.hh"
#include "utils.hh"

Validation::Validation()
{
}

/// Sets up a tornado run
void Validation::tornado_setup(string add_info, Operation mode, ExtFactor ext_factor, string file, const vector <string> &data_sim_lines) 
{
	auto nc = num_core();
	if(nc != 1) run_error("Setup requires only one core");
	
	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	
	auto name = file.substr(0,file.length()-5);
	
	auto tor_dir = name+"_tornado/";
	auto run_dir = tor_dir+"run/";
	auto op_dir = tor_dir+"op/";
	
	auto ncore_inf = 1u;
	
	auto nrun = TORNADO_NUM;
	auto para_run = true;
	
	if(add_info.length() > 0){
		add_info = remove_quote(add_info);
		auto spl = split(add_info,' ');
		
		for(const auto &va: spl){
			if(begin_str(va,"n=")){
				nrun = number(va.substr(2));
				if(nrun == UNSET) alert_input("After 'tornado' the value "+va+" should set a number");
			}
			else{
				if(va == "para") para_run = true;
				else{
					if(va == "series") para_run = false;
					else{
						alert_input("After 'tornado' the value "+va+" is not understood");
					}
				}
			}
		}
	}
	
	auto file_rep = get_file_rep(nrun,run_dir);

	{                                        // Removes data and copies into 
		Model model(mode,ext_factor,true); 
		Mpi mpi(UNSET,model);   
		Output output(model,mpi);    
		Input input(model,file,UNSET,mpi);     // Imports information from input into model

		output.init(input);

		const auto &inf_det = model.inf_details;
		ncore_inf = inf_det.nchain/inf_det.num_per_core;
		
		// Removes any existing data 
		for(auto p = 0u; p < model.nspecies; p++){
			for(const auto &so : model.species[p].source){
				output.delete_command(so.name,so.line_num,false);
			}
		}
		
		output.ensure_directory(tor_dir);
		output.ensure_directory(run_dir);
		output.ensure_directory(op_dir);
		
		{
			ofstream fout(tor_dir+"info.txt");
			fout << "nrun=" << file_rep.size() << endl;
		}
		
		for(auto r = 0u; r < file_rep.size(); r++){
			output.copy(file_rep[r]);
		}
	}
	
	run_sim_data(name,op_dir,para_run,file_rep,ncore_inf,ext_factor,data_sim_lines,file,"tornado-res");
}


/// Gets a list of files for replicates
vector <string> Validation::get_file_rep(unsigned int n, string run_dir) const
{
	vector <string> file_rep;
	for(auto r = 0u; r < n; r++){
		auto fi = run_dir+"Run_"+tstr(r+1)+".bici";
		file_rep.push_back(fi);
	}
	
	return file_rep;
}

	
/// Runs simulates and generates data
void Validation::run_sim_data(string name, string op_dir, bool para_run, vector <string> file_rep, unsigned int ncore_inf, ExtFactor ext_factor, const vector <string> &data_sim_lines, string file, string com) const
{
	percentage_start(RUN_PER);
	for(auto r = 0u; r < file_rep.size(); r++){  
		percentage(r,file_rep.size());
		
		auto fi = file_rep[r]; 
		
		{   // Simulates state
			Model model(SIM,ext_factor,true); 
			Mpi mpi(UNSET,model);   
			Output output(model,mpi,true);    
			Input input(model,fi,r,mpi,true);  
			output.init(input);
			
			Simulate simu(model,output,mpi,true);
			simu.run();
		
			output.end(fi,UNSET);
		}
		
		{ // Generates data
			Model model(DATA_SIM,ext_factor,true); 
			Mpi mpi(UNSET,model);   
			Output output(model,mpi,true);    
			Input input(model,fi,r,mpi,true);  
			output.init(input);
			
			DataSim data_sim(model,output);
			for(auto &ds : data_sim_lines){
				data_sim.run(ds,1);
			}
			
			output.end(fi,UNSET);
		}
	}
	percentage_end();
	
	auto spl = split(name,'/');
	
	auto scr = spl[spl.size()-1]+"_script.sh";             // Output the script for running inference
	ofstream fout(scr);		
	
	for(auto r = 0u; r < file_rep.size(); r++){
		auto fi = file_rep[r]; 
		
		fout << "echo 'Run " << fi << "'" << endl;
		if(para_run) fout << "nohup ";
		if(ncore_inf > 1) fout << "mpirun -n " << ncore_inf << " "; 
		fout << "./bici-para " << fi << " -seed=" << r << " inf";
		if(para_run) fout << " > " << op_dir << "inf_" << r << ".txt 2>&1&";
		fout << endl << endl;
	}
	
	fout << "echo ''" << endl;
	
	//com = replace(com,"\\","\\\\");
	fout << "echo \"WHEN INFERENCE COMPLETE: ./bici-para " << file << " " << com << "\"" << endl; 
	
	cout << endl;
	cout << "MAKE EXECUTABLE: chmod +x " << scr << endl;
	cout << "RUN SCRIPT: ./" << scr << endl; 
}

	
/// Used to order proposal speeds
bool Stat_ord (const Stat &st1, const Stat &st2)                      
{ return (st1.mean < st2.mean); };  


/// Combines results from a tornado run
void Validation::tornado_result(Operation mode, ExtFactor ext_factor, bool no_question, string file)
{
	Model model_base(SIM,ext_factor,no_question);  
	Mpi mpi_base(UNSET,model_base);
	Output op(model_base,mpi_base);      
	Input input(model_base,file,UNSET,mpi_base);
		
	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	auto name = file.substr(0,file.length()-5);
	
	auto tor_dir = name+"_tornado/";
	auto run_dir = tor_dir+"run/";
	auto op_dir = tor_dir+"op/";
	
	auto nrun = TORNADO_NUM;
	
	{
		string info;
		ifstream fin(tor_dir+"info.txt");
		getline (fin,info);
		if(!begin_str(info,"nrun=")) emsg("Problem with loading number of runs");
		nrun = number(info.substr(5));
		if(nrun == UNSET) emsg("Problem with loading number of runs");
	}
	
	auto file_rep = get_file_rep(nrun,run_dir);
		
	vector <Table> param_stats;
	vector <Table> pred_acc_stats;
		
	percentage_start(RUN_PER);
	for(auto r = 0u; r < file_rep.size(); r++){
		percentage(r,file_rep.size());
	
		auto fi = file_rep[r];
		
		Model model(mode,ext_factor,no_question);  
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
	percentage_end();
	
	{
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
			ofstream fout(tor_dir+"ESS.csv");
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
			ofstream fout(tor_dir+"GR.csv");
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
			ofstream fout(tor_dir+tname);
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
					fout << "," << value[i] << "," << stat[i][r].mean << "," << stat[i][r].CImin << "," << stat[i][r].CImax;
				}
				fout << endl;
			}
		}
		
		auto un_name = "uncertainty.csv";
		{	
			ofstream fout(tor_dir+un_name);
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
			ofstream fout(gp+"true_value.csv");
			fout << "y";
			for(auto i = 0u; i < N; i++){
				auto na = add_escape_char(tab0.ele[i][0]);
				fout << ",\"" << na << " value\"";
			}
			fout << endl;
			fout << 0;
			for(auto i = 0u; i < N; i++) fout << "," << value[i];
			fout << endl;
			fout << nrun+1;
			for(auto i = 0u; i < N; i++) fout << "," << value[i];
			fout << endl;
		}
		
		{  // Creates gp scripts
			ofstream fout(gp+"gp_single.txt");
		
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
			ofstream fout(gp+"gp.txt");
		
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
		ofstream fout(tor_dir+tname);
		fout << output_table(tab0);
	}

	cout << endl << "RESULTS: In directory '" << tor_dir << "'." << endl;
	//cout << "PLOT: Gnuplot can be used to plot results. Go to the 'gp' sub-directory and type 'gnuplot gp.txt' to generate graphs and 'ps2pdf gp.ps gp.pdf' to convert to pdf." << endl;
	
	cout << "PLOT: Gnuplot can be used to plot results. Go to the 'gp' sub-directory and type 'gnuplot gp_single.txt' to generate graphs and 'ps2pdf gp_single.ps gp.pdf' to convert to pdf." << endl;
}


/// Sets up a parameter scan
void Validation::scan_setup(string add_info, Operation mode, ExtFactor ext_factor, string file, const vector <string> &data_sim_lines, bool test) 
{
	auto nc = num_core();
	if(nc != 1) run_error("Setup requires only one core");
	
	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	
	auto name = file.substr(0,file.length()-5);
	
	string tor_dir, run_dir, op_dir;

	auto ncore_inf = 1u;
	
	ScanInfo si;
	string na;
	
	vector <string> file_rep;
	{                                        // Removes data and copies into 
		Model model(mode,ext_factor,true); 
		Mpi mpi(UNSET,model);   
		Output output(model,mpi);    
		Input input(model,file,UNSET,mpi);     // Imports information from input into model

		output.init(input);

		si = get_scan_info(add_info,model);

		auto pr = get_param_ref(si.param_name,model);
		if(pr.th == UNSET || pr.index == UNSET){
			run_error("Parameter '"+si.param_name+"' could not be found in the model.");
		}
		
		na = add_escape_char(si.param_name);
		na = replace(na,"\\","");
		tor_dir = name+"_scan_"+na+"/";
		run_dir = tor_dir+"run/";
		op_dir = tor_dir+"op/";
	
		file_rep = get_file_rep(si.value_sim.size(),run_dir);
		
		const auto &inf_det = model.inf_details;
	
		ncore_inf = inf_det.nchain/inf_det.num_per_core;
		
		// Removes any existing data 
		for(auto p = 0u; p < model.nspecies; p++){
			for(const auto &so : model.species[p].source){
				output.delete_command(so.name,so.line_num,false);
			}
		}
		
		output.ensure_directory(tor_dir);
		output.ensure_directory(run_dir);
		output.ensure_directory(op_dir);
		
		for(auto r = 0u; r < file_rep.size(); r++){
			output.change_sim_value(si.param_name,si.value_sim[r]);
			output.copy(file_rep[r]);
		}
	}
	
	{
		ofstream fout(tor_dir+"info.txt");
		fout << add_info << endl;
	}
	
	run_sim_data(name+"_"+na,op_dir,test,file_rep,ncore_inf,ext_factor,data_sim_lines,file,"scan-res:'"+si.param_name+"'");
}


/// Combines results from a tornado run
void Validation::scan_result(string scan_info, Operation mode, ExtFactor ext_factor, bool no_question, string file)
{
	Model model_base(SIM,ext_factor,no_question);  
	Mpi mpi_base(UNSET,model_base);
	Output op(model_base,mpi_base);      
	Input input(model_base,file,UNSET,mpi_base);

	if(!end_str(file,".bici")) run_error("The file does not end in '.bici'");
	auto name = file.substr(0,file.length()-5);
	
	auto spl = split(scan_info,' ');
	auto pname = spl[0];
	
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
	
	string sc_inf;
	ifstream fin(tor_dir+"info.txt");
	getline (fin,sc_inf);
	
	auto si = get_scan_info(sc_inf,model_base);

	auto prop = get_param_prop(si.param_name);
	
	auto file_rep = get_file_rep(si.value_sim.size(),run_dir);
		
	vector <Table> param_stats;
	
	percentage_start(RUN_PER);
	for(auto r = 0u; r < file_rep.size(); r++){
		percentage(r,file_rep.size());
	
		auto fi = file_rep[r];
		
		Model model(mode,ext_factor,no_question);  
		Mpi mpi(UNSET,model);
		Output output(model,mpi,true);           
		Input input(model,fi,UNSET,mpi,true);    
		
		const auto &tab = model.inf_param_stats;
		if(tab.heading.size() == 0) run_error("Parameter statistics could not be found. Have you run the script?");
		
		param_stats.push_back(tab);
	}
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
		ofstream fout(tor_dir+"ESS.csv");
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
		ofstream fout(tor_dir+"GR.csv");
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
	
	auto isel = UNSET;
	auto tname = "result.csv";
	{	
		ofstream fout(tor_dir+tname);
		fout << "Run";
		for(auto i = 0u; i < N; i++){
			auto na = add_escape_char(tab0.ele[i][0]);
			auto prop2 = get_param_prop(na);
			if(prop2.name == prop.name){
				auto d = 0u; while(d < prop.dep.size() && prop.dep[d] == prop2.dep[d]) d++;
				if(d == prop.dep.size()) isel = i;
			}
			fout << ",\"" << na << " value\"" << ",\"" << na << " (post. mean)\"" << ",\"" << na << " (CI min)\"" << ",\"" << na << " (CI max)\"";
		}
		fout << endl;
		if(isel == UNSET) emsg("Cannot find selected parameter");
		
		for(auto r = 0u; r < R; r++){
			fout << (r+1);
			for(auto i = 0u; i < N; i++){
				
				fout << ",";
				if(i == isel) fout << si.value_sim[r];
				else fout << value[i];
				fout << "," << stat[i][r].mean << "," << stat[i][r].CImin << "," << stat[i][r].CImax;
			}
			fout << endl;
		}
	}

	auto gp = tor_dir+"gp/";
	op.ensure_directory(gp);
	
	{  // Creates gp scripts
		ofstream fout(gp+"gp_single.txt");
	
		fout << gp_header();
		if(si.log_on){
			fout << "set logscale x" << endl; 
			fout << "set logscale y" << endl; 
		}
		fout << "set size square" << endl;
		
		auto namex = "True "+add_escape_char(tab0.ele[isel][0]);
		namex = convert_label(namex);
			
		for(auto i = 0u; i < N; i++){
			auto namey = add_escape_char(tab0.ele[i][0]);
			namey = convert_label(namey);
			
			fout << "set xlabel '"+namex+"' font ',30'" << endl;
			fout << "set ylabel '"+namey+"' font ',30'" << endl;
			
			fout << "plot '../" << tname << "' using " << 2+4*isel << ":" << 3+4*i << " ps 1 lw 4 lc '#000000' notitle,\\" << endl;
			fout << "'../" << tname << "' using " << 2+4*isel << ":" << 3+4*i << ":" << 3+4*i+1 << ":" << 3+4*i+2 << " with yerrorbars ls 1 lw 4 notitle,\\" << endl;
			fout << "'../" << tname << "' using " << 2+4*isel << ":" << 2+4*i << " with lines ls 2 notitle" << endl;
		}
	}
	
	{  // Creates gp scripts
		ofstream fout(gp+"gp.txt");
	
		fout << gp_header_multi();
		if(si.log_on){
			fout << "set logscale x" << endl; 
			fout << "set logscale y" << endl; 
		}
		fout << "set size square" << endl;
		
		auto namex = "True "+add_escape_char(tab0.ele[isel][0]);
		namex = convert_label(namex);
			
		for(auto i = 0u; i < N; i++){
			auto namey = add_escape_char(tab0.ele[i][0]);
			namey = convert_label(namey);
			
			fout << "set xlabel '"+namex+"' font ',15'" << endl;
			fout << "set ylabel '"+namey+"' font ',15'" << endl;
			
			fout << "plot '../" << tname << "' using " << 2+4*isel << ":" << 3+4*i << " ps 1 lw 2 lc '#000000' notitle,\\" << endl;
			fout << "'../" << tname << "' using " << 2+4*isel << ":" << 3+4*i << ":" << 3+4*i+1 << ":" << 3+4*i+2 << " with yerrorbars ls 1 lw 2 notitle,\\" << endl;
			fout << "'../" << tname << "' using " << 2+4*isel << ":" << 2+4*i << " with lines ls 2 lw 2 notitle" << endl;
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


/// Gets tha value of a parameter
double Validation::get_param_value(string na, const Model &model) const 
{
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


/// Gets information for a scan info
ScanInfo Validation::get_scan_info(string scan_info, const Model &model) const
{
	auto spl = split(scan_info,' ');
	auto na = spl[0];

	ScanInfo si;
	si.param_name = spl[0];
	
	auto n = 20u;
	auto log_on = false;//true;
	
	auto val = get_param_value(na,model);
	
	double min = val/4, max = val*4;

	for(auto i = 1u; i < spl.size(); i++){
		auto spl2 = split(spl[i],'=');
		if(spl2.size() != 2){
			run_error("Syntax error in scan information '"+scan_info+"'.");
		}
		if(spl2[0] == "range"){
			auto spl3 = split(spl2[1],',');
			
			auto fl = false;
			if(spl3.size() != 2) fl = true;
			else{
				min = number(spl3[0]);
				max = number(spl3[1]);
				if(min == UNSET || max == UNSET) fl = true;
			}
			
			if(fl){
				run_error("In scan information '"+scan_info+"' the value for range '"+spl2[1]+"' is not valid.");
			}
		}
		else{
			if(spl2[0] == "log"){
				if(spl2[1] == "true") log_on = true;
				else{
					if(spl2[1] == "false") log_on = false;
					else{
						run_error("In scan information '"+scan_info+"' the value for log '"+spl2[1]+"' is not valid.");
					}
				}
			}
			else{
				if(spl2[0] == "n"){
					n = number(spl2[1]);
					if(n == UNSET) run_error("In scan information '"+scan_info+"' the value for n '"+spl2[1]+"' is not valid.");
					if(n <= 0) run_error("In scan information '"+scan_info+"' the value for n '"+spl2[1]+"' must be positive.");
					if((unsigned int) n != n) run_error("In scan information '"+scan_info+"' the value for n '"+spl2[1]+"' must be an integer.");
				}
				else{
					run_error("In scan information '"+scan_info+"' the value '"+spl2[0]+"' is not understood.");
				}
			}
		}
	}		

	if(min >= max) run_error("In scan information '"+scan_info+"' the minimum must be smaller than the maximum.");
	
	if(log_on){
		if(min <= 0) run_error("In scan information '"+scan_info+"' the minimum must be positive when log is used.");
	}

	for(auto i = 0u; i <= n; i++){
		if(log_on){
			si.value_sim.push_back(exp(log(min)+(log(max)-log(min))*double(i)/n));	
		}
		else{
			si.value_sim.push_back(min+(max-min)*double(i)/n);
		}
	}
	
	si.log_on = log_on;
	
	return si;
}
