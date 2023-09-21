#ifndef BICI__INPUT_HH
#define BICI__INPUT_HH

#include <string>

using namespace std;

#include "model.hh"

struct Fragment {                              // A fragment of an input line 
	string text;
	unsigned int pos, pos_end;
	unsigned int quote;
};

struct Tag {                                   // Stores information about a single tag in a command line
	string name;
	unsigned int pos;
	unsigned int pos_end;
	string value;
	unsigned int done;
};

struct CommandLine {                           // Stores a command line instruction
	Command command;
	string command_name;
	unsigned int type_pos;
	vector <Tag> tags;
	unsigned int line_num;
};

struct Bracket {                               // Stores bracket information
	string ch;  
	unsigned int i;
};

struct DataTemplate {                    // Used to store information about loading data columns
	DataTemplate(Command cna, vector <string> co){ cname = cna; cols = co;};
	Command cname;                         // The heading for the column
	vector <string> cols;                  // Columns needed for data source
};
	
const vector <DataTemplate> data_template = {
DataTemplate(INIT_POP_PRIOR,{"cl_all","alpha"}),
DataTemplate(INIT_POP,{"cl_all","pop"}),
DataTemplate(ADD_IND,{"ID","t,the time the individuals are added","cl_all"}),
DataTemplate(REMOVE_IND,{"ID","t,the time the individuals are removed"}),
DataTemplate(MOVE_IND,{"ID","t,the time the individuals are moved","to"}),
DataTemplate(COMP_DATA,{"ID","t,the time the individual's state is measured","cl_prob"}),
DataTemplate(TRANS_DATA,{"ID","t,the time the individuals undergo the transition"}),
DataTemplate(TRANS_TIMERANGE_DATA,{"ID","t,the time the individuals undergo the transition","start","end"}),
DataTemplate(SOURCE_DATA,{"ID","t,the time the individuals enter the system"}),
DataTemplate(SOURCE_TIMERANGE_DATA,{"ID","t,the time the individuals enter the system","start","end"}),
DataTemplate(SINK_DATA,{"ID","t,the time the individuals leave the system"}),
DataTemplate(SINK_TIMERANGE_DATA,{"ID","t,the time the individuals leave the system","start","end"}),
DataTemplate(TEST_DATA,{"ID","t,the time the individual's diagnostic test is taken","result"}),
DataTemplate(POP_DATA,{"t,the time the population measurement is taken","filt_obspop"}),
DataTemplate(POP_TRANS_DATA,{"t,the time the individuals undergo the transition","filt_obspoptrans"}),
DataTemplate(SET_TRAPS_DATA,{"setrap_name","t,the time the population measurement is taken","filt_settraps"}),
DataTemplate(IND_TRAPPED_DATA,{"ID","setrap_name"}),
DataTemplate(GENETIC_DATA,{"ID","t,the time the genetic measurement is taken","snp,the first SNP column"}),
};
	
	
class Input                                           // Stores information about the model
{
	public:
		string outputdir;                                 // The output directory 
		vector <string> lines_raw;                        // Stores the raw text lines of the input file
		
		Input(Model &model, string file);
	
	private:
		unsigned int p_current, cl_current;               // Stores the current species and classification 
		
		bool terminate;                                   //
		
		string input_file;                                // The input file
		
		string datadir;                                   // The data directory
		
		vector <string> lines;                            // Stores the text lines of the input file
		
		unsigned int line_num;                            // Stores the line being processed
	
		CommandLine cline_store;                          // Stores the line for the current command being processed
		string tag_find;                                  // The current tag being found
		
		vector <GeoJSON> geo_json;                        // Stores information about loaded geoJSON file
		
		vector <ErrorMess> error_mess;                    // Stores error messages
		
		Model &model;                                     // References model such that it can be updated
		
		vector <CommandLine> extract_command_line(vector <string> lines);
		string add_escape_char(string te);
		CommandLine get_command_tags(string trr, unsigned int line_num);
		CommandLine syntax_error() const;
		void alert(string st);
		void alert_import(string st);  
		void alert_line(string st, unsigned int line);  
		void alert_equation(const EquationInfo &eqi, const string &warn);
		void alert_warning(string st);  
		bool fatal_error() const;
		void process_command(const CommandLine &cline, unsigned int loop);
		string get_tag_value(string st);
		void cannot_find_tag();
		unsigned int option_error(string na, string te, const vector <string> &pos, const vector <unsigned int> &conv);
		bool check_latlng_error(double lat, double lng);
		void output_error_messages(string te) const;
		bool is_number(string num, string tag);
		bool is_positive(string num, string tag);
		bool is_zeroone(string num, string tag);
		Prior convert_text_to_prior(string te);
		void map_ind_effect();
		void check_param_used();
		void create_param_vector();
		void create_pop_ref();
		void raw_process();
		//vector <PriorPos> set_pri_pos(type);
		
		// In 'input_commands.cc' 
		void species_command(unsigned int loop);
		void classification_command(unsigned int loop);
		void set_command();
		void camera_command();
		void compartment_command();
		void transition_command(const Command cname);
		void datadir_command();
		void description_command();
		void label_command();
		void box_command();
		void param_command();
		void derived_command();
		void simulation_command();
		void inference_command();
		void ind_effect_command();
		void fixed_effect_command();
		void import_data_table_command(Command cname);
		void dummy_file_command();
		
		// In 'input_modelupdate.cc'
		void add_species(string name, SpeciesType sp_type);
		void add_classification(unsigned int p, string name, string index, Coord coord);
		void add_compartment(string name, unsigned int p, unsigned int cl, double x, double y, double lat, double lng, bool markov_branch, string erlang_source);
		void add_transition(unsigned int p, unsigned int cl, unsigned int i,unsigned int f, TransType type);
		void determine_branching() const;
		EquationInfo add_equation_info(string _te, EqnType _type, unsigned int _p = UNSET, unsigned int _cl = UNSET);
		void create_equations();
		void add_parent_child(const EquationInfo eqi, unsigned int i, unsigned int th);
		void create_markov_eqn();
		void create_markov_eqn_pop_ref();
		void global_comp_trans_init();
		void create_population_erlang();
		string get_erlang_name(string name, const vector <Compartment> &comp) const;
		void branch_param_group();
		void population_bp_rate_combine();
		void source_equation_comp();
		void calculate_timepoint();	
		void create_spline();
		bool is_dist(const Param &par) const;
		void set_dist(Param &par);
		void create_nm_trans();
		void ind_eff_group_trans_ref();
		void ind_eff_pop_ref();
		void param_affect_likelihood();
		void create_markov_comp_gl();
		void add_pop_affect(unsigned int k, unsigned int po);
		void create_island();
		void create_trg_from_tr();
		
		// In 'input_utils.cc'
		unsigned int find_p(string name) const;
		unsigned int find_cl(unsigned int p, string name) const;
		unsigned int find_cl_index(unsigned int p, string index) const;
		unsigned int find_c(unsigned int p, unsigned int cl, string name) const;
		unsigned int find_tr(unsigned int p, unsigned int cl, string name) const;
		//void clone_compartments(unsigned int p_from, unsigned int cl_from, unsigned int p_to, unsigned int cl_to);
		bool check_claa_error();
		Classification& get_claa();
		bool check_char_allowed(string st, string not_allowed);
		bool check_comp_exist(string name, unsigned int p);
		Table load_table(const string file);
		Table get_subtable(const Table &tab, const vector <string> &col_name);
		unsigned int get_cl_from_comp(string name, unsigned int p) const;
		unsigned int find_string_in(const vector <string> &arr, string val) const;
		bool set_loadcol(Command cname, DataSource &ds);
		void load_obs_model(ObsModel &om);
		string stringify(const vector <string> &arr) const;
		CompRef find_comp_from_name(unsigned int p, string te) const;
		unsigned int import_geojson(string file);
		LatLng boundary_mean_latlng(unsigned int i, string name);
		//ParamProp get_param_prop(string st);
		void set_element(vector <EquationInfo> &value, const vector <Dependency> &dep, const vector <unsigned int> &ind, string te);
		void set_reparam_element(vector <EquationInfo> &value, const vector <Dependency> &dep, const vector <unsigned int> &ind, const EquationInfo &val);
		void set_prior_element(vector <Prior> &prior, const vector <Dependency> &dep, const vector <unsigned int> &ind, Prior pri);
		void clone_class(unsigned int p_to, unsigned int p_from, unsigned int cl_from);
		void print_table(const Table &tab) const;
		void add_to_list(vector <ParamRef> &list, const ParamRef &pr) const;
		unsigned int get_dependency(vector <Dependency> &dep, const ParamProp &pp, const vector <string> &knot_times);

		// In 'input_check.cc'
		void check_initial_pop_error(bool end);
		void data_source_check_error(const DataSource &ds);
		string check_element(const DataSource &ds, unsigned int r, unsigned int c) const;
		Result check_eqn_valid(string te) const;
		void check_comp_structure();
		void check_import_correct();
		void temp_check(unsigned int num);
		string basic_equation_check(string te);
		unsigned int param_end(const string &st, unsigned int i, string &warn) const;
		bool check_brackets_match(string te);
		bool check_chnotallowed(string te);
		unsigned int check_pos_integer(string val, unsigned int def);
		double check_pos_number(string te, unsigned int def);
		double check_zero_one(string te, double def);

		// In 'input_data.cc'
		//void initialise_data();
};
#endif