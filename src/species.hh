#ifndef BICI__SPECIES_HH
#define BICI__SPECIES_HH

#include <vector>

using namespace std;

#include "struct.hh"
#include "equation.hh"
//#include "model.hh"

class Species   	                       // Stores details of a species
{ 	     
	public:
	Species(const vector <double> &timepoint, const Details &details);
	

		string name;                           // The name of the species
		SpeciesType type;                      // Type of model (individual / population-based)

		vector <Classification> cla;           // Stores classification within the species
		unsigned int ncla;                     // Number of classifications

		unsigned int N;                        // The total number of compartments in species

		vector <IndEffect> ind_effect;         // Any individual effects within species
		
		vector <IEgroup> ind_eff_group;        // Groups of individual effects
		
		vector <FixedEffect> fix_effect;       // Any fixed effects within species
		
		vector <unsigned int> comp_mult;       // Used to convert from classification definition to global

		vector <CompGlobal> comp_gl;           // Global compartments

		vector <TransGlobal> tra_gl;           // Global transitions  
		
		vector < vector < vector <unsigned int> > >  trg_from_tr; // Converts from tr to trg [cgl][cl][tr]
		
		vector <NMTrans> nm_trans;             // A list of non-Markovian transitions
		
		vector <MarkovEqn> markov_eqn;         // List Markovian equations
			
		MarkovTree markov_tree;                // Constructs a tree of nodes to sample markov equations
		
		vector <IndData> individual;           // Stores data on individuals
			
		vector <DataSource> source;            // Sources of data
		InitCond init_cond;                    // Stores infomation about the initial conditions
		
		vector <unsigned> obs_eqn;             // Lists observation equations
		
		vector <PopData> pop_data;             // Stores population data
		
		vector <PopTransData> pop_trans_data;  // Stores population transition data
		
		vector < vector <double> > trans_not_allow; // This capture transitions not allowed by transition data
		
		vector < vector < vector <double> > > pop_trans_ref; // References transition
		
		vector <WarnData> warn;                // Stores any errors used to generate data

		unsigned int T;                                     // The maximum time
		
		void create_markov_tree();
		unsigned int find_c(unsigned int cl, string name) const;
		unsigned int update_c_comp(unsigned int c, unsigned int cl, unsigned int c_comp) const;
		//void print_event_data(string name, const Individual &ind);
		void print_obs_data(string name, const vector <ObsData> &obs) const;
		unsigned int get_ti(double t) const;
	
	// In 'species_data.cc':
	public: 
		void initialise_data();
		void order_data_events();
		
 	private:
		void init_pop_data(const DataSource &so);
		void add_ind_data(const DataSource &so);
		void remove_ind_data(const DataSource &so);
		void move_ind_data(const DataSource &so);
		void trans_data(const DataSource &so);
		void comp_data(const DataSource &so);
		void test_data(const DataSource &so);
		void population_data(const DataSource &so);
		void population_trans_data(const DataSource &so);
		vector <bool> global_convert(const Filter &filt) const;
		vector <bool> trans_global_convert(unsigned int cl, const vector <bool> &trans_filt, const Filter &filt);
	
		unsigned int find_individual(string name);
		void X_vector_order();
		
		
		void generate_data() const;
		void alert_source(string st, const DataSource &so, unsigned int c = UNSET, unsigned int r = UNSET);
	public:
	const vector <double> &timepoint;
	const Details &details;
};
#endif
