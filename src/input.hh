#pragma once

#include <string>

using namespace std;

#include "model.hh"
#include "hash.hh"
#include "mpi.hh"

struct TransDef {                          // Transition reference
	TransType type;
	string rate;
	string mean;
	string cv;
	string shape;
	string scale;
	string time;
	bool set;
};

struct Fragment {                          // A fragment of an input line 
	string text;                             // Text in fragment
	unsigned int pos, pos_end;               // Position of fragment
	unsigned int quote;                      
};

struct FileStore {                         // Used to store a fi;e
	string name;                             // Filename
	char sep;                                // Seperator
	vector <string> lines;                   // Lines in file
};

struct Tag {                               // Stores information about a single tag in a command line
	string name;
	unsigned int pos;
	unsigned int pos_end;
	string value;
	unsigned int done;
};

struct CommandLine {                       // Stores a command line instruction
	Command command;
	string command_name;
	unsigned int type_pos;
	vector <Tag> tags;
	unsigned int all_row;
	unsigned int line_num;
};

struct DataTemplate {                      // Used to store information about loading data columns
	DataTemplate(Command cna, vector <string> co){ cname = cna; cols = co;};
	Command cname;                           // The heading for the column
	vector <string> cols;                    // Columns needed for data source
};
	

const vector <DataTemplate> data_template = {
DataTemplate(INIT_POP,{"cl_all","init_pop"}),
DataTemplate(ADD_POP,{"t,the time the individuals are added","cl_all","add_pop"}),
DataTemplate(REMOVE_POP,{"t,the time the individuals are removed","cl_all","rem_pop"}),
DataTemplate(ADD_IND,{"ID","t,the time the individuals are added","cl_all_prob"}),
DataTemplate(REMOVE_IND,{"ID","t,the time the individuals are removed"}),
DataTemplate(MOVE_IND,{"ID","t,the time the individuals are moved","to"}),
DataTemplate(COMP_DATA,{"ID","t,the time the individual's state is measured","cl_prob"}),
DataTemplate(TRANS_DATA,{"ID","t,the time the individuals undergo the transition"}),
DataTemplate(TRANS_TIMERANGE_DATA,{"ID","t,the time the individuals undergo the transition","start","end"}),
DataTemplate(TEST_DATA,{"ID","t,the time the individual's diagnostic test is taken","result"}),
DataTemplate(POP_DATA,{"t,the time the population measurement is taken","filt_obspop"}),
DataTemplate(POP_TRANS_DATA,{"tstart","tend","filt_obspoptrans"}),
DataTemplate(GENETIC_DATA,{"ID","t,the time the genetic measurement is taken","snp,the first SNP column"}),
};
	
class Input                                // Stores information about the model
{
	public:
		string datadir;                        // The data directory 
		vector <string> lines_raw;             // Stores the raw text lines of the input file
		
		Input(Model &model, string file, unsigned int seed, Mpi &mpi);
	
	private:
		unsigned int p_current, cl_current;    // Stores the current species and classification 
		
		bool terminate;                        // Set if algorithm terminated early
		
		string input_file;                     // The input file
	
		vector <FileStore> files;              // Stores files from the input file
	
		vector <string> lines;                 // Stores the text lines of the input file
		
		unsigned int line_num;                 // Stores the line being processed
	
		CommandLine cline_store;               // Stores the line for the current command being processed
		string tag_find;                       // The current tag being found
		unsigned int all_row;                  // The row when comp_all or trans_all are loaded
		
		vector <GeoJSON> geo_json;             // Stores information about loaded geoJSON file
		
		vector <ErrorMess> error_mess;         // Stores error messages
		
		Hash hash_eqn;                         // Stores a hash take for equations
	
		Model &model;                          // References model such that it can be updated
		Mpi &mpi;                              // Stores information about mpi
		
		vector <CommandLine> extract_command_line(vector <string> lines);
		void load_data_files(vector <CommandLine> &command_line);
		string add_escape_char(string te);
		void convert_folder(string &data_dir) const;
		CommandLine get_command_tags(string trr, unsigned int line_num);
		CommandLine syntax_error() const;
		void alert(string st);
		void alert_import(string st, bool fatal=false);  
		void alert_line(string st, unsigned int line);  
		void alert_equation(const EquationInfo &eqi, const string &warn);
		void alert_warning(string st);  
		bool fatal_error() const;
		void process_command(const CommandLine &cline, unsigned int loop);
		string get_tag_value(string st);
		string get_tag_val(string st, vector <Tag> &tags);
		void cannot_find_tag();
		unsigned int option_error(string na, string te, const vector <string> &pos, const vector <unsigned int> &conv);
		bool check_latlng_error(double lat, double lng);
		void output_error_messages(string te) const;
		bool is_number(string num, string tag);
		bool is_positive(string num, string tag);
		bool is_zeroone(string num, string tag);
		void map_ind_effect();
		void check_param_used();
		void create_param_vector();
		void create_pop_ref();
		void simplify_equations();
		
		// In 'input_commands.cc' 
		void species_command(unsigned int loop);
		void classification_command(unsigned int loop);
		void param_mult_command();
		void set_command();
		void camera_command();
		vector < vector <Tag> > get_tags_list(string file);
		void check_tags_used(unsigned int r, const vector <Tag> &tags);
		void compartment_command();
		void compartment_all_command();
		void compartment_command2(vector <Tag> &tags);
		void transition_command();
		void transition_command2(vector <Tag> &tags);
		void transition_all_command();
		void datadir_command();
		void description_command();
		void label_command();
		void box_command();
		void param_command();
		void load_reparam_eqn(string te, Param &par);
		void derived_command();
		void simulation_command();
		void inference_command();
		void post_sim_command();
		void ind_effect_command();
		void fixed_effect_command();
		void import_data_table_command(Command cname);
		void map_command();
		void inf_state_command();
		void dummy_file_command();
		
		// In 'input_modelupdate.cc'
		void add_species(string name, SpeciesType sp_type, bool trans_tree);
		void add_classification(unsigned int p, string name, string index, Coord coord);
		void add_compartment(string name, unsigned int p, unsigned int cl, double x, double y, double lat, double lng, bool markov_branch, CompInfected infected, string erlang_source);
		void add_transition(unsigned int p, unsigned int cl, unsigned int i,unsigned int f, TransType type);
		void determine_branching() const;
		void add_obs_model_eqn(Species &sp);
		void create_equations();
		void add_eq_ref(EquationInfo &eqi);	
		vector <bool> set_eqn_zero(const vector <EquationInfo> &eq_info);
		void add_parent_child(const EquationInfo eqi, unsigned int i, unsigned int th);
		void create_markov_eqn();
		void create_markov_eqn_pop_ref();
		void global_comp_trans_init();
		void set_tr_leave_markov();
		void set_multi_trans_cl();
		void set_tr_connected();
		void set_precalc_nm_rate();
		void check_markov_or_nm() const;
		void check_nm_pop() const;
		void create_population_erlang();
		void source_equation_comp();
		void calculate_timepoint();	
		void create_spline();
		bool is_dist(const Param &par) const;
		void create_nm_trans();
		void set_tra_ie_fe_affect();
		void ind_fix_eff_group_trans_ref();
		void ind_fix_eff_pop_ref();
		void param_affect_likelihood();
		void create_markov_comp_gl();
		void omr_add_to_map(vector <unsigned int> &omr, unsigned int i, vector < vector <bool> > &map, const vector < vector <bool> > &spline_map, const Species &sp) const;
		void add_nm_trans_affect(unsigned int p, unsigned int i, unsigned int eq, AffectType type, const vector < vector <bool> > &spline_map);
		void add_obs_trans_eqn(unsigned int p, unsigned int e, const vector < vector <bool> > &spline_map);
		void add_to_map(unsigned int eq, unsigned int i, vector < vector <bool> > &map, const vector < vector <bool> > &spline_map) const;
		void add_map_list(const vector < vector <bool> > &map, unsigned int p, AffectType type) const;
		void create_island();
		void create_trg_from_tr();
		void source_rate_divide() const;
		void setup_obs_trans_is_one();
		void setup_obs_trans();		
		void bp_create_unset();
		void markov_bp_convert();
		void exp_nm_convert();
		void setup_trans_infection();
		void set_trans_tree();
		void set_hash_all_ind();
		void add_genetic_data();
		void set_contains_source();
		void create_trig_ev_ref();
		bool trig_ev_exist(TrigEventType type, double t, const vector <TrigEventRef> &trig_event) const;
		void create_cl_trig_ev_ref();
		void set_joint_param_event();
		void set_local_ind_init();
		void print_tr_swap(unsigned int cl, const vector < vector <unsigned int> > &tr_swap_ref, const vector < vector < vector <TrSwap> > > &tr_swap_list, const Species &sp) const;
		bool add_tr_swap(unsigned int ci, unsigned int cf, const vector <TrSwap> &vec, vector < vector <unsigned int> > &tr_swap_leave, vector < vector <unsigned int> > &tr_swap_ref, vector < vector < vector <TrSwap> > > &tr_swap_list);
		void set_comp_global_convert();
		void set_tra_global_time_vari();
		void set_ind_variation();
		void set_init_c_set();
		void set_cgl_begin_nm();
		void set_cgl_tr_source_sink();   
		void set_tr_enter();
		void set_inf_cause();
		void set_ppc_resample();
		void set_comp_period();
		void set_comp_terminal();
		void set_sink_exist();
		void set_eqn_ind_eff_exist();
		void set_param_parent_child();
		void set_param_use();
		void set_omega_fl();
	
		// In 'input_utils.cc'
		unsigned int find_p(string name) const;
		unsigned int find_cl(unsigned int p, string name) const;
		unsigned int find_cl_index(unsigned int p, string index) const;
		unsigned int find_c(unsigned int p, unsigned int cl, string name) const;
		bool check_claa_error();
		Classification& get_claa();
		bool check_char_allowed(string st, string not_allowed);
		bool check_comp_exist(string name, unsigned int p);
		Table load_table(const string file);
		Table get_subtable(const Table &tab, const vector <string> &col_name);
		unsigned int get_cl_from_comp(string name, unsigned int p) const;
		unsigned int get_cl_from_trans(string name, unsigned int p) const;
		unsigned int find_string_in(const vector <string> &arr, string val) const;
		bool set_loadcol(Command cname, DataSource &ds);
		void set_SNP_columns(const Table &tab, DataSource &ds);
		void set_genetic_matrix_columns(const Table &tab, DataSource &ds);
		void load_obs_model(ObsModel &om);
		string stringify(const vector <string> &arr) const;
		CompRef find_comp_from_name(unsigned int p, string te) const;
		unsigned int import_geojson(string file);
		LatLng boundary_mean_latlng(unsigned int i, string name);
		void set_element(vector <EquationInfo> &value, const vector <Dependency> &dep, const vector <unsigned int> &ind, string te);
		void set_reparam_element(vector <EquationInfo> &value, const vector <Dependency> &dep, const vector <unsigned int> &ind, const EquationInfo &val);
		void set_prior_element(vector <Prior> &prior, const vector <Dependency> &dep, const vector <unsigned int> &ind, Prior pri);
		void clone_class(unsigned int p_to, unsigned int p_from, unsigned int cl_from);
		void print_table(const Table &tab) const;
		void add_to_list(vector <ParamRef> &list, const ParamRef &pr) const;
		unsigned int get_dependency(vector <Dependency> &dep, const ParamProp &pp, const vector <string> &knot_times);
		EquationInfo he(EquationInfo eqn_inf);
		bool is_file(string te) const;
		vector <unsigned int> find_index(unsigned int i, const vector <Dependency> &depend) const;
		TransDef extract_trans_def(string value) const;
		string get_prop(string value, string prop, string end) const;
		string in_file_text(string te) const;
		void alert_sample(string warn, unsigned int num);
		void read_state_sample(const vector <string> &lines, const vector <string> &ind_key);
		unsigned int get_param_value(vector < vector <double> > &param_value, unsigned int i, const vector <string> &lines, string warn);
		void load_param_value(const ParamProp &pp, string valu, Param &par, string desc);
		void set_spline(string knot_times_str, string smooth, vector <string> &knot_times, bool use_inf_time, Param &par);
		unsigned int get_seed();
		//unsigned int get_tag_integer(string tag, unsigned int def=UNSET);
		
		// In 'input_check.cc'
		void check_initial_pop_error(bool end);
		void data_source_check_error(const DataSource &ds);
		string check_element(const DataSource &ds, unsigned int r, unsigned int c);
		Result check_eqn_valid(string te) const;
		void check_comp_structure();
		void check_import_correct();
		void temp_check(unsigned int num);
		unsigned int check_pos_integer(string val, unsigned int def=UNSET);
		double check_pos_number(string te, unsigned int def);
		double check_zero_one(string te, double def);
};

