/// General purpose utility functions for all classes to use

#ifndef BICI__UTILS_HH
#define BICI__UTILS_HH

using namespace std;

#include "struct.hh"

void emsg(const string &msg);
string trim(string te);
string remove_prime(string te);
string remove_quote(string te);
vector<string> split(const string &s, char delimiter);
string toLower(string st);
string toUpper(string st);
unsigned int find_in(const vector <unsigned int> &vec, const unsigned int val);
unsigned int find_in(const vector <string> &vec, const string val);
unsigned int find_in(const vector <string> &vec, const string val, unsigned int guess);
unsigned int find_in(string st1, string st2);
unsigned int find_in(const vector <bool> &vec, const bool val);
double number(string st);
double ran();
void sran(const int seed);
double normal_sample(const double mean, const double sd);
double normal_probability(const double x, const double mean, const double sd);
double lognormal_sample(const double mean, const double cv);
double lognormal_probability(const double x, const double mean, const double cv); 
double weibull_sample(const double shape, const double scale);
double weibull_probability(const double x, const double shape, const double scale);
double gamma_sample(const double mean, const double cv);
double gamma_probability(const double x, const double mean, const double cv);
double beta_sample(const double alpha, const double beta);
double beta_probability(const double x, const double alpha, const double beta);
unsigned int poisson_sample(const double lam);
double poisson_probability(const int i, const double lam);
unsigned int bernoulli_sample(const double p);
double bernoulli_probability(unsigned int x, const double p);
double exp_rate_sample(const double rate);
double exp_rate_probability(const double x, const double rate); 
double exp_mean_sample(const double mean);
double exp_mean_probability(const double x, const double mean);
double gamma_alpha_sample(const double alpha);
double gamma_alpha_probability(const double x, const double alpha);
unsigned int binomial_sample(double p, unsigned int n);
double binomial_probability(unsigned int num, double p, unsigned int n);

void test_distribution();
double rectify(double val);

bool equal_vec(const vector <string> &vec1, const vector <string> &vec2);
bool equal_vec(const vector <unsigned int> &vec1, const vector <unsigned int> &vec2);
string stringify(const vector <EquationInfo> &value, const vector <Dependency> &dep);
bool includes(string st, string ch);
DepInfo get_dep_info(const string &te, unsigned int i, const string &not_allowed);
bool str_eq(const string &te, unsigned int i, const string &te2);
SwapResult swap_index(string &te, const vector <DepConv> &dep_conv);
ParamProp get_param_prop(string st);
unsigned int add_to_vec(vector <unsigned int> &vec, unsigned int val);
void print(string name, const vector <double> &vec);
void print(string name, const vector <unsigned int> &vec);
void printsm(string name, const vector <double> &vec);
void print(string name, const vector < vector <double> > &vec);
void print(string name, const Table &tab);
bool dif(double a, double b);
bool dif(const vector < vector <double> > &a, const vector < vector <double> > &b);

string replace(string st, string st1, string st2);
void remove_cr(string &st);
bool is_percent(string val);
string get_param_name_with_dep(const Param &par, const vector <Dependency> &dep, unsigned int index);
void param_vec_add_affect(vector <AffectLike> &vec, const AffectLike &al);
double sum(const vector <double> &vec);
double sum(const vector < vector <double> > &mat);
string cpu_percent(long time, long total_time);
string tstr(double value);
void progress(double val, double val2);
#endif
