// Stores functions for sampling from different distributions and deals with error messages

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <stdexcept>
#include "stdlib.h"
#include "math.h"
#include <sys/stat.h>
#include <cstring>
#include <signal.h>
#include <math.h>
//#include <numbers> 
//#include "mpi.hh"

#include <algorithm>
//#include <cctype>
//#include <iomanip>

using namespace std;

#include "const.hh"
#include "utils.hh"

default_random_engine generator;

/// Displays an error message
void emsg(const string &msg)
{
	//if(emsg_throws) throw(std::runtime_error(msg));
	if(com_op == true){
		cout << "<<ERROR>>" << endl <<  msg << endl << "<<END>>" << endl;
	}
	else{
		cout << "\033[31m" << "ERROR:" << endl;
		cout << "\033[0m";
		
		cout << msg;
		if(msg.length() > 0 && msg.substr(msg.length()-1,1) != ".") cout << ".";
		cout << endl;
		
		
		//if(debug == true) raise(SIGABRT);
	}
	
	exit (EXIT_FAILURE);
}


/// Trims white space from a string
string trim(string te)
{
	auto i = 0u; while(i < te.length() && te.substr(i,1) == " ") i++;
	if(i > 0) te = te.substr(i);

	if(te.length() > 0){
		i = te.length()-1; while(i > 0 && te.substr(i,1) == " ") i--;
		if(i < te.length()-1) te = te.substr(0,i+1);
	}
	
	return te;
}	


/// Removes any primes
string remove_prime(string te)
{
	if(te == "") return "";
	int i = te.length()-1;

	while(i >= 0 && te.substr(i,1) == "'") i--;
	return te.substr(0,i+1);
}


/// Removes quotation marks from a string
string remove_quote(string te)
{
	te = trim(te);
	if(te.length() >= 2 && te.substr(0,1) == "\"" && te.substr(te.length()-1,1) == "\"") te = te.substr(1,te.length()-2);
	if(te.length() >= 2 && te.substr(0,1) == "'" && te.substr(te.length()-1,1) == "'") te = te.substr(1,te.length()-2);
	return te;
}


/// Split up a string at a specified delimiter
vector<string> split(const string &s, char delimiter)                                                          
{                              
  vector<string> splits;                       
 
  bool quoteon = false;
	auto j = 0u;
	for(auto i = 0u; i < s.length(); i++){
		if(s.substr(i,1) == "\""){
			if(quoteon == false) quoteon = true; else quoteon = false;
		}
		
		if(s.at(i) == delimiter && quoteon == false){
			splits.push_back(s.substr(j,i-j)); j = i+1; 
		}
	}
	splits.push_back(s.substr(j,s.length()-j));
	for(auto &spl : splits) spl = trim(spl);
	
	return splits;                                           
}


/// Coverts a string to lower case
string toLower(string st)
{	
	transform(st.begin(), st.end(), st.begin(), ::tolower);
	return st;
}


/// Coverts a string to upper case
string toUpper(string st)
{	
	transform(st.begin(), st.end(), st.begin(), ::toupper);
	return st;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <unsigned int> &vec, const unsigned int val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <string> &vec, const string val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <string> &vec, const string val, unsigned int guess)
{
	if(guess < vec.size() && vec[guess] == val) return guess;
	
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Finds the index of a value in a vector
unsigned int find_in(string st1, string st2)
{
	auto i = st1.find(st2);
	if(i == string::npos) return UNSET;
	return i;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <bool> &vec, const bool val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == val) return i;
	}
	return UNSET;
}


/// Checks a string has a certain set of characters
bool allow_string(const string st, const string ok_char)
{
	for(auto i = 0u; i < st.length(); i++){
		auto j = 0u; while(j < ok_char.length() && st.substr(i,1) != ok_char.substr(j,1)) j++;
		if(j == ok_char.length()) return false;
	}
	return true;
}


/// Converts a string to a number
double number(string st)
{
	double val = atof(st.c_str());
	if(val == 0 && trim(st) != "0" && trim(st).substr(0,3) != "0.0") return UNSET;
	return val;

/*
st = "H";
	stringstream ss(st);
	double val;
	ss >> val;

	cout << st << " " << val << "  val\n";

	return val;
*/
/*
	st = trim(st);
	if(allow_string(st,"-0123456789.eE") == false) return UNSET;
	else{
		char* endptr;
		double val = strtod(&st[0],&endptr);
		ptrdiff_t j = endptr-&st[0];
		if(j != (int)st.length()) return UNSET;
		return val;
	}
*/

}



/// Sets the seed for the random number generator
void sran(const int seed)
{
	generator.seed(seed);
	srand(seed);
	//mt = std::mt19937(seed);
}


/// Draws a random number between 0 and 1
double ran()
{
	return double(0.999999999*rand())/RAND_MAX;
	//std::uniform_real_distribution<> dist(0.0000000001,0.999999999);
	//return dist(mt);
}


/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(const double mean, const double sd) 
{
	normal_distribution<double> distribution(mean, sd);
	return distribution(generator);
}


/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double sd)
{
	auto var = sd*sd;
  if(var <= 0) return -LARGE;
  return -0.5*log(2*M_PI*var) - (x-mean)*(x-mean)/(2*var);
}


/// Draws a log-normally distributed number with mean mu and standard deviation sd
double lognormal_sample(const double mean, const double cv) 
{
	if(cv <= 0) emsg("Utils");
	if(mean <= 0) emsg("Utils");
	
	auto var = log(1+cv*cv);   // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	normal_distribution<double> distribution(mu, sqrt(var));
	
	return exp(distribution(generator));
}


// Probability of log-normally distribution
double lognormal_probability(const double x, const double mean, const double cv) 
{
	if(mean <= 0 || cv <= 0 || x <= 0) return -LARGE;
	
	auto var = log(1+cv*cv);   // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	auto logx = log(x);
	return  -0.5*log(2*M_PI*var*x*x) - (logx-mu)*(logx-mu)/(2*var);
}


/// Draws a Weibull distributed number with shape and scale parameters
double weibull_sample(const double shape, const double scale) 
{
	if(shape <= 0 || scale <= 0) emsg("Utils");
	
	weibull_distribution<double> distribution(shape, scale);
	return distribution(generator);
}


/// Probability of weibull sample
double weibull_probability(const double x, const double shape, const double scale) 
{
	if(shape <= 0 || scale <= 0) return -LARGE;
	
	return log(shape/scale)+(shape-1)*log(x/scale) - pow(x/scale,shape);
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
double gamma_sample(const double mean, const double cv)
{
	gamma_distribution<double> distribution(1.0/(cv*cv),mean*cv*cv);
	return distribution(generator);
}


/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double mean, const double cv)
{
	if(mean < TINY || cv < TINY) return -LARGE;
	
	// Scale parameter is 1/b in several sources
	auto shape = 1.0/(cv*cv);
	auto b = shape/mean;
  if(x <= 0 || shape <= 0 || b <= 0) return -LARGE;
  return (shape-1)*log(x) - b*x + shape*log(b) - lgamma(shape);
}


/// Draws a sample from the beta distribution x^(a-1)*(1-x)^(b-1)
double beta_sample(const double alpha, const double beta)
{
	if(alpha < TINY || beta < TINY) emsg("Beta prob");
	
	gamma_distribution<double> distributionX(alpha,1), distributionY(beta,1);
	auto x = distributionX(generator);
	auto y = distributionY(generator);

	return x/(x+y);
}


/// The log of the probability from the beta distribution
double beta_probability(const double x, const double alpha, const double beta)
{
	if(x < TINY || x > 1-TINY) return -LARGE;
	if(alpha < TINY || beta < TINY) return -LARGE;
	
  return (alpha-1)*log(x) + (beta-1)*log(1-x) + lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta);
}


/// Generates a sample from the Poisson distribution
unsigned int poisson_sample(const double lam)
{
	if(lam > LARGE) emsg("Utils");
  poisson_distribution<int> distribution(lam);
	return distribution(generator);
}


/// The log probability of the poisson distribution
double poisson_probability(const int i, const double lam)
{
	if(lam < 0) return -LARGE;
	
  if(lam == 0){
    if(i == 0) return 0;
    else return -LARGE;
  }
  else{
    if(i < 0) return -LARGE;
    else{
			return i*log(lam)-lam - lgamma(i+1);
		}
  }
}


/// Generates a sample from the Bernoulli distribution
unsigned int bernoulli_sample(const double p)
{
	if(p < 0) emsg("Utils");
		
	if(ran() < p) return 1;
	return 0;
}


/// The probability of a Bernoulli sample
double bernoulli_probability(unsigned int x, const double p)
{
	if(p <= 0 || p >= 1) return -LARGE;
	if(x == 0) return log(1-p);
	return log(p);
}


/// Samples from the exponential distribution with specified rate
double exp_rate_sample(const double rate)
{
	return -log(ran())/rate;
}


/// Probability of exponential distribution
double exp_rate_probability(const double x, const double rate)
{
	return log(rate) - rate*x;
}

/// Samples from the exponential distribution with specified mean
double exp_mean_sample(const double mean)
{
	return -log(ran())*mean;
}


/// Probability of exponential distribution
double exp_mean_probability(const double x, const double mean)
{
	return -log(mean) - x/mean;
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-x)
double gamma_alpha_sample(const double alpha)
{
	gamma_distribution<double> distribution(alpha,1);
	return distribution(generator);
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-x)
double gamma_alpha_probability(const double x, const double alpha)
{
  if(x <= 0 || alpha <= 0) return -LARGE;
  return (alpha-1)*log(x) - x - lgamma(alpha);
}


/// A sample from the binomial distribution
unsigned int binomial_sample(double p, unsigned int n)
{
	if(n == 0) return 0;
		
	binomial_distribution<int> distribution(n,p);
	return distribution(generator);
}


/// A sample from the binomial distribution
double binomial_probability(unsigned int num, double p, unsigned int n)
{
	if(num > n || p < 0 || p > 1) emsg("Problem");
	
	if(p == 0){
		if(num == 0) return 0;
		else return -LARGE;
	}
	else{
		return lgamma(n+1) - lgamma(num+1) - lgamma(n-num+1) + num*log(p) + (n-num)*log(1-p);
	}
}


/// Ensures a quantity is positive or zero
double rectify(double val)
{
	if(val < 0) return 0;
	return val;
}


/// Compares two strings to see if they are equal
bool equal_vec(const vector <string> &vec1, const vector <string> &vec2)
{
	if(vec1.size() != vec2.size()) return false;
	
	for(auto i = 0u; i < vec1.size(); i++){
		if(vec1[i] != vec2[i]) return false;
	}
	
	return true;
}



/// Compares two strings to see if they are equal
bool equal_vec(const vector <unsigned int> &vec1, const vector <unsigned int> &vec2)
{
	if(vec1.size() != vec2.size()) return false;
	
	for(auto i = 0u; i < vec1.size(); i++){
		if(vec1[i] != vec2[i]) return false;
	}
	
	return true;
}


/// Converts a value tensor into a string
string stringify(const vector <EquationInfo> &value, const vector <Dependency> &dep)
{
	string st = "";
	for(auto i = 0u; i < value.size(); i++){
		auto num_begin = 0u, num_end = 0u; 
		for(auto j = 0u; j < dep.size(); j++){
			auto si = dep[j].mult*dep[j].list.size();
			if(i%si == 0) num_begin++;
			if(i%si == si-1) num_end++;
		}

		for(auto k = 0u; k < num_begin; k++) st += "[";
		st += value[i].te;
		if(num_end == 0) st += ",";
		else st += "]";
	}

	return st;
}


/// Determinets if a string contains a character
bool includes(string st, string ch)
{
	auto loc = st.find(ch);
	if(loc != string::npos) return true;
	
	return false;
}


/// Gets the dependency from an initial point on the string
DepInfo get_dep_info(const string &te, unsigned int i, const string &not_allowed)
{
	DepInfo di; di.result = SUCCESS;

	string dep;

	auto ist = i;
	if(te.substr(i,1) == "("){	
		while(i < te.length() && te.substr(i,1) != ")") i++;
		if(i == te.length()){ di.result = FAIL; return di;}

		dep = te.substr(ist+1,i-(ist+1));
	}
	else{
		while(i < te.length() && !includes(not_allowed,te.substr(i,1)) && te.substr(i,1) != "[") i++;
		dep = te.substr(ist,i-ist);
	}

	if(dep == ""){
		//di.result = FAIL;
	}
	else{
		auto spl = split(dep,',');

		for(auto &de : spl){
			de = trim(de);
			if(de == ""){ di.result = FAIL; return di;}
		}

		di.spl = spl;
	}
	di.iend = i;

	return di;
}


/// Determines if te2 is at the location i in te
bool str_eq(const string &te, unsigned int i, const string &te2)
{
	if(i+te2.length() > te.length()) return false;
	if(te.substr(i,te2.length()) == te2) return true;
	return false;
}


/// Swaps indices with comparmtent values
SwapResult swap_index(string &te, const vector <DepConv> &dep_conv)
{
	SwapResult res; res.warn = ""; res.done.resize(dep_conv.size(),false);

	auto i = 0u;
	while(i < te.length()){
		while(i < te.length() && te.substr(i,1) != "%" && te.substr(i,1) != "{") i++;
		if(i < te.length()){
			auto istore = i; 

			string type = "";
			if(te.substr(i,1) == "{"){  // Changes in a population
				type = "pop";
				i++;
				auto ist = i;
				while(i < te.length() && te.substr(i,1) != ":" &&  te.substr(i,1) != "}") i++;
				if(i < te.length() && te.substr(i,1) == ":") i++;
				else i = ist;
			}
			else{                       // Changes in a parameter
				while(i < te.length() && te.substr(i,1) != "_" && te.substr(i,1) != "$") i++;
		
				if(i < te.length() && te.substr(i,1) == "_"){
					type = "param";
					i++;
				}
			}
		
			if(type != ""){
				auto ist = i;
				auto di = get_dep_info(te,i,"$}[<(");

				if(di.result == FAIL){
					res.warn = "An error occured";
					if(type == "pop"){
						res.warn = "Population {} has an error";
					}
					if(type == "param"){
						res.warn = "Parameter '"+te.substr(istore+1,i-istore-2)+"' has a misspecified dependency.";
					}
					return res;
				}
		
				string dep_new = "";
				for(auto k = 0u; k < di.spl.size(); k++){
					auto de = di.spl[k];
					for(auto kk = 0u; kk < dep_conv.size(); kk++){
						const auto &dp = dep_conv[kk];
						if(de == dp.before){ de = string(dp.after); res.done[kk] = true;}
					}
				
					if(dep_new != "") dep_new += ",";
					dep_new += de;
				}

				te = te.substr(0,ist)+dep_new+te.substr(di.iend);
				i = ist+dep_new.length()+2;
			}
		}
	}

	return res;
}


/// Gets the properties of a parameter from the string
ParamProp get_param_prop(string st)
{
	st = trim(st);
	auto time_dep = false, age_dep = false;
	
	if(st.length() > 3){
		auto end = st.substr(st.length()-3,3);
		if(end == "(t)") time_dep = true;
		if(end == "(a)") age_dep = true;
		if(end == "(t)" || end == "(a)") st = st.substr(0,st.length()-3);
	}
	
	auto j = 0u;
	while(j < st.length() && st.substr(j,1) != "_") j++;

	auto name = st.substr(0,j);

	vector <string> dep;
	auto prime = false;
	
	if(st.substr(j,1) == "_"){
		j++;
		auto jst = j; 
		while(j < st.length() && st.substr(j,1) != "[") j++;
		auto sub = st.substr(jst,j-jst);
		dep = split(sub,',');
	}
					
	if(time_dep == true) dep.push_back("t");
	if(age_dep == true) dep.push_back("a");
	
	auto dep_with_prime = dep;
		
	for(auto &de : dep){
		auto de_st = de;
		de = remove_prime(de);
		if(de != de_st) prime = true;
	}
		
	ParamProp pp; 
	pp.name = name;
	pp.dep = dep;
	pp.dep_with_prime = dep_with_prime;
	pp.time_dep = time_dep;
	pp.age_dep = age_dep;
	pp.prime = prime;
	
	return pp;
}


/// Adds a value to a vector (if it doesn't already exist)
unsigned int add_to_vec(vector <unsigned int> &vec, unsigned int val)
{
	auto i = find_in(vec,val);
	if(i == UNSET){ i = vec.size(); vec.push_back(val);}
	return i;
}


/// Prints a vector
void print(string name, const vector <double> &vec)
{
	cout << "Vector " << name << ":" << endl;
	for(auto i = 0u; i < vec.size(); i++){
		if(i != 0) cout << ", ";
		cout << vec[i];
	}		
	cout << endl;
}


/// Prints a vector
void print(string name, const vector <unsigned int> &vec)
{
	cout << "Vector " << name << ":" << endl;
	for(auto i = 0u; i < vec.size(); i++){
		if(i != 0) cout << ", ";
		cout << vec[i];
	}		
	cout << endl;
}

/// Prints a vector on a single line
void printsm(string name, const vector <double> &vec)
{
	cout << name << ": ";
	for(auto i = 0u; i < vec.size(); i++){
		if(i != 0) cout << ", ";
		cout << vec[i];
	}		
	cout << endl;
}


/// Prints a matrix
void print(string name, const vector < vector <double> > &mat)
{
	cout << "Matrix " << name << ":" << endl;
	for(auto j = 0u; j < mat.size(); j++){
		for(auto i = 0u; i < mat[j].size(); i++){
			if(i != 0) cout << ", ";
			cout << mat[j][i];
		}		
		cout << endl;
	}
}


/// Prints a table
void print(string name, const Table &tab)
{
	cout << "Table " << name << ":" << endl;
	for(auto i = 0u; i < tab.ncol; i++){
		if(i != 0) cout << ",";
		cout << tab.heading[i];
	}
	cout << endl;
	
	for(auto j = 0u; j < tab.nrow; j++){
		for(auto i = 0u; i < tab.ncol; i++){
			if(i != 0) cout << ",";
			cout << tab.ele[j][i];
		}
		cout << endl;
	}
	cout << endl;
}


/// Calculates if two numbers is different (subject to numerical noise) 
bool dif(double a, double b)
{
	auto dif = a-b; if(dif < 0) dif = -dif;
	//if(dif > TINY) return true;
#ifdef WINDOWS 
	if(dif > 10){
#else
	if(dif > 0.00001){
#endif
		//cout << a << " "<< b << " Different values\n";
		return true;
	}
	return false;
}


/// Calculates if two numbers is different (subject to numerical noise) 
bool dif(const vector < vector <double> > &a, const vector < vector <double> > &b)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(a[j].size() != b[j].size()) return true;
		for(auto i = 0u; i < a[j].size(); i++){
			if(dif(a[j][i],b[j][i])) return true;
		}
	}
	return false;
}



/// Replaces one substring with another
string replace(string st, string st1, string st2)
{
	if(st.length() < st1.length()) return st;
	
	auto i = 0u;
	while(i <= st.length()-st1.length()){
		if(st.substr(i,st1.length()) == st1) st = st.substr(0,i)+st2+st.substr(i+st1.length());
		else i++;
	}
	return st;
}


/// Removes any \r at the end of a line
void remove_cr(string &st)
{
	string end = "\ufeff"; 
	auto len = end.length();
	if(st.length() >= len){
		if(st.substr(0,len) == end) st = st.substr(len);
	}
			
	if(st.length() > 0){
		if(st.substr(st.length()-1) == "\r") st = st.substr(0,st.length()-1);
	}
}


/// Detemines if the value is a percentage
bool is_percent(string val)
{
	if(val.length() == 0) return false;
	if(val.substr(val.length()-1,1) == "%") return true;
	return false;
}


/// Gets the name of a parameter along with a specified dependency
string get_param_name_with_dep(const Param &par, const vector <Dependency> &dep, unsigned int index)
{
	auto name = par.name;
	if(dep.size() > 0){
		name += "_";
		for(auto d = 0u; d < dep.size(); d++){
			if(d != 0) name += ",";
			name += dep[d].list[(index/dep[d].mult)%dep[d].list.size()];
		}	
	}
	
	return name;
}


/// Converts from type to string
string get_affect_name(AffectType type)
{
	switch(type){
	case LIKE_OBS_AFFECT: return "LIKE_OBS_AFFECT";
	case MARKOV_POP_AFFECT: return "MARKOV_POP_AFFECT";
	case LIKE_IE_AFFECT: return "LIKE_IE_AFFECT";
	case EXP_IE_AFFECT: return "EXP_IE_AFFECT";
	case OMEGA_AFFECT: return "OMEGA_AFFECT";
	case POP_AFFECT: return "POP_AFFECT";
	case SPLINE_PRIOR_AFFECT: return "SPLINE_PRIOR_AFFECT";
	case PRIOR_AFFECT: return "PRIOR_AFFECT";
	case DIST_AFFECT: return "DIST_AFFECT";
	case SPLINE_AFFECT: return "SPLINE_AFFECT";
	case EXP_FE_AFFECT: return "EXP_FE_AFFECT";
	case INDFAC_INT_AFFECT: return "INDFAC_INT_AFFECT";
	case DIV_VALUE_AFFECT: return "DIV_VALUE_AFFECT";
	case MARKOV_LIKE_AFFECT: return "MARKOV_LIKE_AFFECT";
	case NM_TRANS_AFFECT: return "NM_TRANS_AFFECT";
	}
	emsg("affect type not known");
	return "error";
}


/// Adds a likelihood affect to a parameter in param_vec 
void param_vec_add_affect(vector <AffectLike> &vec, const AffectLike &al)
{
	switch(al.type){
	case SPLINE_AFFECT: case DIV_VALUE_AFFECT: case MARKOV_LIKE_AFFECT: case POP_AFFECT:
	case MARKOV_POP_AFFECT:
		if(al.map.size() == 0) emsg("map should be set: "+get_affect_name(al.type)); 
		break;
		
	default:
		if(al.map.size() != 0) emsg("map should not be set: "+get_affect_name(al.type));
		break;
	}
	
	auto i = 0u;
	while(i < vec.size() && !(vec[i].type == al.type && vec[i].num == al.num && vec[i].num2 == al.num2)) i++;

	if(i == vec.size())	vec.push_back(al);
	else{
		if(vec[i].map.size() != al.map.size()) emsg("map sizes do not agree");
		
		for(auto t = 0u; t < vec[i].map.size(); t++){
			if(al.map[t] == true) vec[i].map[t] = true;
		}			
	}
	
	vec[i].list.clear();
	for(auto t = 0u; t < vec[i].map.size(); t++){
		if(vec[i].map[t] == true) vec[i].list.push_back(t);
	}
}


/// Sums up the values in a vector
double sum(const vector <double> &vec)
{
	auto sum = 0.0; for(auto val : vec) sum += val;
	return sum;
}


/// Sums up the values in a matrix
double sum(const vector < vector <double> > &mat)
{
	auto su = 0.0; for(const auto &vec : mat) su += sum(vec);
	return su;
}


/// Converts from clock to percent usage
string cpu_percent(long time, long total_time)
{ 
	stringstream ss;
	ss << (unsigned int)((100.0*time)/total_time) << "%";
	return ss.str();
}


/// Converts a value to a string
string tstr(double value)
{
	stringstream ss; ss << value; 
	return ss.str();
}


/// Outputs progess (when doing command line output)
void progress(double val, double val2)
{
	cout << "<PROGRESS>" << (unsigned int)(100*(val+1)/val2) << endl;
	cout.flush();
}


