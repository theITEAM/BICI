// Stores functions for sampling from different distributions and deals with error messages

// The relative speeds for different sampling/probability functions
// lognormal sample - 22%
// lognormal probability - 12%
// lognormal upper probability - 16%
// weibull sample - 15%
// weibull probability - 13%
// weibull upper probability - 7%
// gamma sample - 18%
// gamma probability - 10%
// gamma upper probability - 100%
// exp_rate sample - 6%
// exp_rate probability - 3%
// exp_rate upper probability - 0%
// exp_mean sample - 6%
// exp_mean probability - 7%
// exp_mean upper probability - 0%

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
#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

#include "const.hh"
#include "utils.hh"

default_random_engine generator;
std::uniform_real_distribution<> dist(0,1);
std::random_device rd; 
std::mt19937 gen(rd());

vector <double> log_sum;
vector <double> log_integer;

/// Displays an error message
void emsg(const string &msg)
{
	if(false) throw(std::runtime_error(msg));
	
	display_error(msg);
	if(false) raise(SIGABRT);
	
	exit (EXIT_FAILURE);
}


/// Displays an error message
void display_error(const string &msg)
{
	if(com_op == true){
		cout << "<<ERROR>>" << endl <<  msg << endl << "<<END>>" << endl;
	}
	else{
		cout << "\033[31m" << "ERROR:" << endl;
		cout << "\033[0m";
		
		cout << msg;
		if(msg.length() > 0 && msg.substr(msg.length()-1,1) != ".") cout << ".";
		cout << endl;
	}
}


/// Display a warning message
void display_warning(const string &msg)
{
	if(com_op == true){
	}
	else{
		cout << "\033[35m";
		cout << "WARNING: ";
		cout << "\033[0m";
		cout << msg;
		if(msg.length() > 0 && msg.substr(msg.length()-1,1) != ".") cout << ".";
		cout << endl;
	}
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


/// Split up a string accounting for brackets and quotation marks
vector<string> split_with_bracket(const string &s, char delimiter)
{                              
  vector<string> splits;                       
 
  bool quoteon = false;
	auto num_brac = 0;
	
	auto j = 0u;
	for(auto i = 0u; i < s.length(); i++){
		auto letter = s.substr(i,1);
		
		if(letter == "\""){
			if(quoteon == false) quoteon = true; else quoteon = false;
		}
	
		if(letter == "(") num_brac++;
		if(letter == ")") num_brac--;
		
		if(s.at(i) == delimiter && num_brac == 0 && quoteon == false){
			splits.push_back(s.substr(j,i-j)); j = i+1; 
		}
	}
	splits.push_back(s.substr(j,s.length()-j));
	for(auto &spl : splits) spl = trim(spl);
	
	if(false){
		for(auto va : splits){
			cout << va << " split" << endl;
		}
	}
	
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
unsigned int find_in(const vector <double> &vec, const double val)
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
	if(allow_string(st,"-0123456789.e") == false) return UNSET;
	
	double val = atof(st.c_str());
	if(val == 0 && trim(st) != "0" && trim(st).substr(0,3) != "0.0") return UNSET;
	return val;
}


/// Gets an integer from a string
unsigned int integer(string st)
{
	auto val = number(st);
	if(val == UNSET) return UNSET;
	
	auto i = (unsigned int)(val);
	if(val != i) return UNSET;
	
	return i;
}


/// Sets the seed for the random number generator
void set_seed(const int chain, const Details &details, unsigned int seed_tag)
{
	auto seed = details.seed;
	if(seed == UNSET) seed = SEED_DEFAULT;

	// Seed from tag loading BICI overides seed from file
	if(seed_tag != UNSET) seed = seed_tag;
	
	seed += 10000*chain;

	generator.seed(seed);
	srand(seed);
}


/// Draws a random number between 0 and 1
double ran()
{
	double val;
	do{
		//val = double(rand())/RAND_MAX;
		val = dist(gen);
	}while(val == 0 || val == 1);

	return val;
}


/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(const double mean, const double sd) 
{
	if(sd < 0) emsg("Normal distribution cannot have negative SD");
	normal_distribution<double> distribution(mean, sd);
	return distribution(generator);
}


/// Returns an integer samping distribution with a normal envelope
int normal_int_sample(double si)
{
	return int(normal_sample(0,si)+0.5+LARGE_INT)-LARGE_INT;
}


/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double sd)
{
	auto var = sd*sd;
  if(var < TINY) return -LARGE;
  return -0.5*log(2*M_PI*var) - (x-mean)*(x-mean)/(2*var);
}


/// Draws a log-normally distributed number with mean mu and standard deviation sd
double lognormal_sample(const double mean, const double cv) 
{
	if(cv < 0) emsg("Lognormal cannot have negative CV");
	if(mean < TINY) emsg("Lognormal must have positive mean");
	
	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	normal_distribution<double> distribution(mu, sqrt(var));
	
	return exp(distribution(generator));
}


/// Probability of log-normally distribution
double lognormal_probability(const double x, const double mean, const double cv) 
{
	if(mean < TINY || cv < 0 || x < TINY) return -LARGE;
	
	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	auto logx = log(x);
	return  -0.5*log(2*M_PI*var*x*x) - (logx-mu)*(logx-mu)/(2*var);
}


/// The lognormal probability of being xmin or above
double lognormal_upper_probability(const double xmin, const double mean, const double cv)
{
	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	auto val = 1-0.5*erfc(-(log(xmin)-mu)/(sqrt(2*var)));
	if(val < TINY) val = TINY;
	return log(val);
}


/// The lognormal probability of being xmin or above (with no log taken)
double lognormal_upper_probability_no_log(const double xmin, const double mean, const double cv)
{
	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	auto val = 1-0.5*erfc(-(log(xmin)-mu)/(sqrt(2*var)));
	if(val < TINY) val = TINY;
	return val;
}


/// Draws a Weibull distributed number with shape and scale parameters
double weibull_sample(const double scale, const double shape) 
{
	if(scale <= TINY) emsg("Weibull must have posivie scale");
	if(shape <= TINY) emsg("Weibull must have positive shape");
	
	weibull_distribution<double> distribution(shape, scale);
	return distribution(generator);
}


/// Probability of weibull sample
double weibull_probability(const double x, const double scale, const double shape) 
{
	if(shape < TINY || scale < TINY) return -LARGE;
	
	return log(shape/scale)+(shape-1)*log(x/scale) - pow(x/scale,shape);
}


/// Probability of weibull sample beign xmin or above
double weibull_upper_probability(const double xmin, const double scale, const double shape) 
{
	
	return -pow((xmin/scale),shape);
}


/// Probability of weibull sample beign xmin or above (with no log taken)
double weibull_upper_probability_no_log(const double xmin, const double scale, const double shape) 
{
	
	return exp(-pow((xmin/scale),shape));
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
double gamma_sample(const double mean, const double cv)
{
	if(mean < TINY) emsg("Gamma mean must be positive");
	if(cv < TINY) emsg("Gamma CV must be positive");
	gamma_distribution<double> distribution(1.0/(cv*cv),mean*cv*cv);
	return distribution(generator);
}


/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double mean, const double cv)
{
	if(mean < TINY || cv < TINY) return -LARGE;
	
	auto shape = 1.0/(cv*cv);
	auto b = shape/mean;
  if(x <= 0 || shape <= 0 || b <= 0) return -LARGE;
  return (shape-1)*log(x) - b*x + shape*log(b) - lgamma(shape);
}


/// The gamma probability of being x or above
double gamma_upper_probability(const double xmin, const double mean, const double cv)
{
	auto shape = 1.0/(cv*cv);
	auto theta = mean/shape;
	return log(boost::math::gamma_q(shape,xmin/theta));
}


/// The gamma probability of being x or above (with no log taken)
double gamma_upper_probability_no_log(const double xmin, const double mean, const double cv)
{
	auto shape = 1.0/(cv*cv);
	auto theta = mean/shape;
	return boost::math::gamma_q(shape,xmin/theta);
}


/// Draws a sample from the beta distribution x^(a-1)*(1-x)^(b-1)
double beta_sample(const double alpha, const double beta)
{
	if(alpha < TINY) emsg("For the Beta distribution alpha must be positive");
	if(beta < TINY) emsg("For the Beta distribution beta must be positive");
	
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
	if(lam > LARGE) emsg("For the Poisson distribution the mean is too large");
	if(lam < 0) emsg("Poisson mean cannot be negative");
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


/// The log probability of the negative binomial distribution
double neg_binomial_probability(const int k, const double mean, const double p)
{	
	auto r = round_int(mean*p/(1-p));
	if(r < 1) r = 1;
	return choose(k+r-1,k) + k*log(1-p) + r*log(p);
}


/// Generates a sample from the Bernoulli distribution
unsigned int bernoulli_sample(const double p)
{
	if(p < 0) emsg("Bernoulli mean cannot be smaller than zero");
	if(p > 1) emsg("Bernoulli mean cannot be larger than one");
		
	if(ran() < p) return 1;
	return 0;
}


/// The probability of a Bernoulli sample
double bernoulli_probability(unsigned int x, const double p)
{
	
	if(x == 0){
		if(p > 1-TINY) return -LARGE;
		return log(1-p);
	}
	
	if(p < TINY) return -LARGE;
	return log(p);
}


/// Samples from the exponential distribution with specified rate
double exp_rate_sample(const double rate)
{
	if(rate < TINY) emsg("Exponential rate must be positive");
	return -log(ran())/rate;
}


/// Probability of exponential distribution
double exp_rate_probability(const double x, const double rate)
{
	if(rate < TINY) return -LARGE;
	 
	return log(rate) - rate*x;
}


/// Probability of exponential distribution above xmin
double exp_rate_upper_probability(const double xmin, const double rate)
{
	return -rate*xmin;
}


/// Probability of exponential distribution above xmin (with no log taken)
double exp_rate_upper_probability_no_log(const double xmin, const double rate)
{
	return exp(-rate*xmin);
}


/// Samples from the exponential distribution with specified mean
double exp_mean_sample(const double mean)
{
	if(mean < 0) emsg("The exponential must be positive");
	return -log(ran())*mean;
}


/// Probability of exponential distribution
double exp_mean_probability(const double x, const double mean)
{
	if(x < 0 || mean < TINY) return -LARGE;
	return -log(mean) - x/mean;
}


/// Probability of exponential distribution above xmin
double exp_mean_upper_probability(const double xmin, const double mean)
{
	return -xmin/mean;
}


/// Probability of exponential distribution above xmin (with no log taken)
double exp_mean_upper_probability_no_log(const double xmin, const double mean)
{
	return exp(-xmin/mean);
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-x)
double gamma_alpha_sample(const double alpha)
{
	if(alpha < TINY) emsg("Alpha must be positive");
	gamma_distribution<double> distribution(alpha,1);
	return distribution(generator);
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-x)
double gamma_alpha_probability(const double x, const double alpha)
{
  if(x <= 0 || alpha < TINY) return -LARGE;
  return (alpha-1)*log(x) - x - lgamma(alpha);
}


/// A sample from the binomial distribution
unsigned int binomial_sample(double p, unsigned int n)
{
	if(n == 0) return 0;
	
	if(p < -TINY) emsg("For the binomial distribution the probability must be postive");
	if(p > OVER_ONE){
		emsg("For the binomial distribution the probability must be one or less");
	}
	
	if(p < TINY) return 0;
	if(p > ALMOST_ONE) return n;
	
	binomial_distribution<int> distribution(n,p);
	return distribution(generator);
}


/// A sample from the binomial distribution
double binomial_probability(unsigned int num, double p, unsigned int n)
{	
	if(p == 0){
		if(num == 0) return 0;
		else return -LARGE;
	}
	else{
		if(p == 1){
			if(num == n) return 0;
			else return -LARGE;
		}
		else{
			if(num > n || p < 0 || p > 1) emsg("Problem");

			return lgamma(n+1) - lgamma(num+1) - lgamma(n-num+1) + num*log(p) + (n-num)*log(1-p);
		}
	}
}


/// Probability of period sample
double period_sample(const double time) 
{
	if(time < TINY) emsg("Period should be positive");
	return time;
}


/// Probability of period sample (Uses a very sharp spike prior)
double period_probability(const double x, const double time) 
{
	if(time < TINY) return -LARGE;
	auto dd = x - time;
	if(dd > 0) return -dd*LARGE;
	else return dd*LARGE;
}


/// Probability of period sample (Uses a very sharp spike prior)
double period_upper_probability(const double x, const double time) 
{
	if(x < time) return 0; 
	return -LARGE;
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
		while(i < te.length() && !includes(not_allowed,te.substr(i,1)) && !str_eq(te,i,sigma) && te.substr(i,1) != "[") i++;
		dep = te.substr(ist,i-ist);
	}

	if(dep != ""){
		auto spl = split(dep,',');

		for(auto &de : spl){
			di.ipos.push_back(ist);
			ist += de.length() + 1;
			
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


/// Swaps indices with comparmtent values (this is now depreciate and only used for checking)
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
			else{                               // Changes in a parameter
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
				i = ist+dep_new.length()+1;
			}
		}
	}

	return res;
}


/// Constructs a swap template
SwapResult swap_template(string te, const vector <DepConv> &dep_conv)
{
	SwapResult res; res.warn = ""; res.done.resize(dep_conv.size(),false);

	auto icopy = 0u;
	
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
			else{                               // Changes in a parameter
				while(i < te.length() && te.substr(i,1) != "_" && te.substr(i,1) != "$") i++;
		
				if(i < te.length() && te.substr(i,1) == "_"){
					type = "param";
					i++;
				}
			}
		
			if(type != ""){
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
		
				for(auto k = 0u; k < di.spl.size(); k++){
					auto de = di.spl[k];
					for(auto kk = 0u; kk < dep_conv.size(); kk++){
						const auto &dp = dep_conv[kk];
						if(de == dp.before){
							auto ipo = di.ipos[k];
							
							SwapTemp swt; 
							swt.te = te.substr(icopy,ipo-icopy);
							swt.num = kk;
							icopy = ipo+de.length();
							res.swap_temp.push_back(swt);
							res.done[kk] = true;
						}
					}
				}
				i = di.iend;
			}
		}
	}
	
	if(icopy < te.length()){
		SwapTemp swt; 
		swt.te = te.substr(icopy);
		swt.num = UNSET;
		res.swap_temp.push_back(swt);
	}
	
	if(false){
		cout << te << " text" << endl;
		for(auto swt : res.swap_temp){
			cout << swt.te << " " << swt.num << endl;
		}
		emsg("Split");
	}
	
	return res;
}


/// Swaps index using a template
string swap_index_temp(const vector <DepConv> &dep_conv, const SwapResult &temp)
{
	string te = "";
	for(const auto &swt : temp.swap_temp){
		te += swt.te;
		if(swt.num != UNSET) te += dep_conv[swt.num].after;
	}
	
	return te;
}


/// Gets the properties of a parameter from the string
ParamProp get_param_prop(string st)
{
	st = trim(st);
	auto time_dep = false;
	
	if(st.length() > 3){
		auto end = st.substr(st.length()-3,3);
		if(end == "(t)"){ time_dep = true;  st = st.substr(0,st.length()-3);}
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


/// Adds a value to a vector (if it doesn't already exist)
unsigned int add_to_vec(vector <ParamRef> &vec, unsigned int th, unsigned int index)
{
	auto i = 0u; 
	while(i < vec.size() && !(vec[i].th == th && vec[i].index == index)) i++;
	if(i == vec.size()){
		ParamRef pref; pref.th = th; pref.index = index; 
		vec.push_back(pref);
	}	
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


/// Prints a matrix
void print(string name, const vector < vector <unsigned int> > &mat)
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
	if(std::isnan(a)) return true;
	if(std::isnan(b)) return true;
	
	auto dif = a-b; if(dif < 0) dif = -dif;

	if(dif > 0.00001){
		if(a < 0) a = -a;
		if(b < 0) b = -b;
		
		auto frac = dif/(a+b);
		if(frac > 0.00000001) return true;
	}
	return false;
}


/// Calculates if two numbers is different (subject to numerical noise) 
bool dif(const vector <double> &a, const vector <double> &b)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(dif(a[j],b[j])) return true;
	}
	return false;
}


/// returns which is different
unsigned int which_dif(const vector <double> &a, const vector <double> &b)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(dif(a[j],b[j])) return j;
	}
	return UNSET;
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


/// Calculates if two numbers is different (subject to numerical noise) 
vector <unsigned int> which_dif(const vector < vector <double> > &a, const vector < vector <double> > &b)
{
	vector <unsigned int> vec;
	
	if(a.size() != b.size()) return vec;
	for(auto j = 0u; j < a.size(); j++){
		if(a[j].size() != b[j].size()) return vec;
		for(auto i = 0u; i < a[j].size(); i++){
			if(dif(a[j][i],b[j][i])){
				vec.push_back(j); vec.push_back(i); 
				return vec;
			}
		}
	}
	return vec;
}


/// Calculates if two unsigned ints are different 
bool dif(unsigned int a, unsigned int b)
{
	if(a != b) return true;
	return false;
}


/// Calculates if two vectors are different 
bool dif(const vector <unsigned int> &a, const vector <unsigned int> &b)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(dif(a[j],b[j])) return true;
	}
	return false;
}


/// Calculates if two numbers is different 
bool dif(const vector < vector <unsigned int> > &a, const vector < vector <unsigned int> > &b)
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
	case OBS_EQN_AFFECT: return "OBS_EQN_AFFECT";
	case LIKE_UNOBS_TRANS_AFFECT: return "LIKE_UNOBS_TRANS_AFFECT";
	case POP_DATA_CGL_TGL_AFFECT: return "POP_DATA_CGL_TGL_AFFECT";
	case LIKE_OBS_IND_AFFECT: return "LIKE_OBS_IND_AFFECT";
	case LIKE_OBS_POP_AFFECT: return "LIKE_OBS_POP_AFFECT";
	case LIKE_OBS_POP_TRANS_AFFECT: return "LIKE_OBS_POP_TRANS_AFFECT";
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
	case DIV_VALUE_FAST_AFFECT: return "DIV_VALUE_FAST_AFFECT";
	case DIV_VALUE_LINEAR_AFFECT: return "DIV_VALUE_LINEAR_AFFECT";
	case MARKOV_LIKE_AFFECT: return "MARKOV_LIKE_AFFECT";
	case NM_TRANS_AFFECT: return "NM_TRANS_AFFECT";
	case NM_TRANS_BP_AFFECT: return "NM_TRANS_BP_AFFECT";
	case NM_TRANS_INCOMP_AFFECT: return "NM_TRANS_INCOMP_AFFECT";
	case LIKE_INIT_COND_AFFECT: return "LIKE_INIT_COND_AFFECT";
	case PRIOR_INIT_COND_AFFECT: return "PRIOR_INIT_COND_AFFECT";
	case LIKE_GENETIC_PROCESS_AFFECT: return "LIKE_GENETIC_PROCESS_AFFECT";
	case GENETIC_VALUE_AFFECT: return "GENETIC_VALUE_AFFECT";
	case LIKE_GENETIC_OBS_AFFECT: return "LIKE_GENETIC_OBS_AFFECT";
	case IIF_W_AFFECT: return "IIF_W_AFFECT";
	case POPNUM_IND_W_AFFECT: return "POPNUM_IND_W_AFFECT";
	case AFFECT_MAX: break;
	}
	emsg("affect type not known");
	return "error";
}


/// Adds a likelihood affect to a parameter in param_vec 
void param_vec_add_affect(vector <AffectLike> &vec, const AffectLike &al)
{
	switch(al.type){
	case SPLINE_AFFECT: case DIV_VALUE_AFFECT: case DIV_VALUE_FAST_AFFECT: 
	case DIV_VALUE_LINEAR_AFFECT: 
	case MARKOV_LIKE_AFFECT: case POP_AFFECT:
	case MARKOV_POP_AFFECT:
	case LIKE_OBS_IND_AFFECT: case LIKE_OBS_POP_AFFECT:
	case LIKE_OBS_POP_TRANS_AFFECT: case OBS_EQN_AFFECT:
	case LIKE_UNOBS_TRANS_AFFECT:
	case NM_TRANS_AFFECT:	case NM_TRANS_BP_AFFECT: case NM_TRANS_INCOMP_AFFECT:
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


/// Sums up the values in a vector
unsigned int sum(const vector <unsigned int> &vec)
{
	auto sum = 0u; for(auto val : vec) sum += val;
	return sum;
}


/// Multiplies elements of matrix and adds
double sum_mult(const vector < vector <double> > &M1, const vector < vector <unsigned int> > &M2)
{
	auto sum = 0.0;
	for(auto j = 0u; j < M1.size(); j++){
		const auto &M1j = M1[j];
		const auto &M2j = M2[j];
		for(auto i = 0u; i < M1j.size(); i++){
			sum += M1j[i]*M2j[i];
		}
	}
	
	return sum;
}


/// Finds the minimum value for an array
double min(const vector <double> &val)
{
	auto ma = LARGE;
	for(auto v : val){
		if(v < ma) ma = v;
	}
	return ma;
}


/// Finds the minimum value for an array
unsigned int min(const vector <unsigned int> &val)
{
	unsigned int ma = LARGE;
	for(auto v : val){
		if(v < ma) ma = v;
	}
	return ma;
}


/// Finds the maximum value for an array
double max(const vector <double> &val)
{
	auto ma = -LARGE;
	for(auto v : val){
		if(v > ma) ma = v;
	}
	return ma;
}

/// Finds the maximum value for an array
unsigned int max(const vector <unsigned int> &val)
{
	auto ma = 0u;
	for(auto v : val){
		if(v > ma) ma = v;
	}
	return ma;
}


/// Rounds a number to the nearest integer 
int round_int(double val)
{
	return int(val+LARGE_INT+0.5)-LARGE_INT;
}


/// Rounds a vector to the nearest integer 
vector <int> round_int(const vector <double> &val)
{
	vector <int> vec(val.size());
	for(auto i = 0u; i < val.size(); i++) vec[i] = round_int(val[i]);
	return vec;
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
	if(com_op && op()){
		cout << "<PROGRESS>" << (unsigned int)(100*(val+1)/val2) << endl;
		cout.flush();
	}
}


/// Determines if output
bool op()
{
#ifdef USE_MPI
	int num;
	MPI_Comm_rank(MPI_COMM_WORLD,&num);
	if(num == 0) return true;
	return false;
#else		
	return true;
#endif
}


/// Determines  number of cores running on
unsigned int num_core()
{
#ifdef USE_MPI
	int num;
	MPI_Comm_size(MPI_COMM_WORLD,&num); 
	return num;
#else		
	return 1;
#endif
}


/// Constructs a sequention vector, e.g. 0,1,2,3,...N-1
vector <unsigned int> seq_vec(unsigned int N)
{
	vector <unsigned int> vec;
	for(auto i = 0u; i < N; i++) vec.push_back(i);
	
	return vec;
}


/// Constructs a vector with all one
vector <bool> true_vec(unsigned int N)
{
	vector <bool> vec(N,true);
	return vec;
}


/// Converts a text string to a prior specification
Prior convert_text_to_prior(string te, unsigned int line_num)
{
	Prior pri;
	
	te = trim(te);
	
	auto spl = split(te,'(');
	if(spl.size() != 2){ pri.error = true; return pri;}
	
	if(spl[1].size() == 0){ pri.error = "Prior has a syntax error"; return pri;}
	if(spl[1].substr(spl[1].length()-1,1) != ")"){ pri.error = "Prior has a syntax error"; return pri;}
	
	auto type = toLower(spl[0]);
	if(type == "dir") type = "dirichlet";
	if(type == "bern") type = "bernoulli";
	
	auto fl = false;
	if(type == "uniform"){ pri.type = UNIFORM_PR; fl = true;}
	if(type == "exp"){ pri.type = EXP_PR; fl = true;}
	if(type == "normal"){ pri.type = NORMAL_PR; fl = true;}
	if(type == "gamma"){ pri.type = GAMMA_PR; fl = true;}
	if(type == "log-normal"){ pri.type = LOG_NORMAL_PR; fl = true;}
	if(type == "beta"){ pri.type = BETA_PR; fl = true;}
	if(type == "bernoulli"){ pri.type = BERNOULLI_PR; fl = true;}
	if(type == "fix"){ pri.type = FIX_PR; fl = true;}
	if(type == "flat"){ pri.type = FLAT_PR; fl = true;}
	if(type == "dirichlet"){ pri.type = DIRICHLET_PR; fl = true;}
	
	if(fl == false){
		pri.error = "Distribution '"+type+"' not recognised"; return pri;
	}
	
	auto bra = spl[1].substr(0,spl[1].length()-1);
	auto spl2 = split(bra,',');

	switch(pri.type){
	case FIX_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
	
	case UNIFORM_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
	
	case EXP_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
		
	case NORMAL_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
		
	case GAMMA_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
		
	case LOG_NORMAL_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;
		
	case BETA_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
		}
		break;

	case BERNOULLI_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
		
	case FLAT_PR:
		break;
		
	case DIRICHLET_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
		}
		break;
	}

	for(auto i = 0u; i < spl2.size(); i++){
		pri.dist_param.push_back(add_equation_info(spl2[i],DIST,line_num));
	}
	
	return pri;
}
	
	
/// Converts a prior to a string
string get_prior_string(Prior prior)
{
	switch(prior.type){
	case FIX_PR:
		return "fix("+prior.dist_param[0].te_raw+")";
		
	case UNIFORM_PR:
		return "uniform("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";
		
	case EXP_PR:
		return "exp("+prior.dist_param[0].te_raw+")";
		
	case NORMAL_PR:
		return "normal("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";
		
	case GAMMA_PR:
		return "gamma("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";
		
	case LOG_NORMAL_PR:
		return "log-normal("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";
		
	case BETA_PR:
		return "beta("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";
	
	case BERNOULLI_PR:
		return "bernoulli("+prior.dist_param[0].te_raw+")";
	
	case FLAT_PR:
		return "flat";
	
	case DIRICHLET_PR:
		return "dir("+prior.dist_param[0].te_raw+")";
	}
	
	return "Problem with prior";
}


/// Intialises the logsum (used to calculate factorials quickly)
void init_log_sum()
{
	log_sum.push_back(0);
	log_integer.push_back(-LARGE);
}


/// Calculate log(n)
double log_int(unsigned int n)
{
	if(n >= log_integer.size()){
		for(auto i = log_integer.size(); i < 2*n; i++){
			log_integer.push_back(log(i));
		}
	}
	
	return log_integer[n];
}


/// Calculate log(n!)
double factorial(unsigned int n)
{
	if(n >= log_sum.size()){
		auto ist = log_sum.size();
		auto sum = log_sum[ist-1];
		
		for(auto i = ist; i < 2*n; i++){
			sum += log(i);
			log_sum.push_back(sum);
		}
	}
	
	return log_sum[n];
}


/// Combinitorics log of the choose function
double choose(unsigned int N, unsigned int r)
{
	return factorial(N) - factorial(N-r) - factorial(r);
}


/// Samples from a Dirichlet distribution
vector <double> dirichlet_sample(const vector <double> &alpha)
{
	auto N = alpha.size();
	
	vector <double> vec(N);
	auto sum = 0.0;
	for(auto i = 0u; i < N; i++){
		auto al = alpha[i];
		if(al == ALPHA_ZERO) vec[i] = 0;
		else{
			vec[i] = gamma_alpha_sample(al);
			sum += vec[i];
		}
	}

	for(auto i = 0u; i < N; i++) vec[i] /= sum;
				
	return vec;
}


/// The probability of a Dirichlet distribution (without constant factor
double dirichlet_probability(const vector <double> &frac, const vector <double> &alpha)
{
	auto N = alpha.size();
	
	auto sum = 0.0;
	for(auto i = 0u; i < N; i++){
		auto al = alpha[i];
		if(al == ALPHA_ZERO) al = 0;
		
		sum += (al-1)*frac[i];
	}
	
	return sum;
}


/// Used to order fractions
bool Frac_ord (FracSort fs1, FracSort fs2)                      
{ return (fs1.frac > fs2.frac); };  


/// Samples from a multinomial probability distribution
vector <unsigned int>  multinomial_sample(unsigned int N_total, const vector <double> &frac)
{
	auto N = frac.size();
	
	vector <unsigned int> x(N,0);
	
	if(false){ // Direct way
		vector <double> frac_sum(N);
		auto sum = 0.0;
		for(auto i = 0u; i < N; i++){
			sum += frac[i];
			frac_sum[i] = sum;
		}
		if(dif(sum,1)) emsg("sum not one");
		
		for(auto j = 0u; j < N_total; j++){
			auto z = ran(); 
			auto i = 0u; while(i < N && z > frac_sum[i]) i++;
			if(i == N) emsg(" Samp problem");
			x[i]++;
		}
	}
	else{  // Faster way
		// Sorts so highest number first
		vector <FracSort> frac_sort;
		for(auto i = 0u; i < N; i++){
			FracSort fs; fs.i = i; fs.frac = frac[i];
			frac_sort.push_back(fs);
		}
		sort(frac_sort.begin(),frac_sort.end(),Frac_ord); 
		
		auto S = N_total;
		auto r = 1.0;
		for(auto k = 0u; k < N-1; k++){
			auto i = frac_sort[k].i;
			auto n = 0.0;
			if(r > 0) n = binomial_sample(ALMOST_ONE*frac[i]/r,S);
			x[i] = n;
			S -= n;
			r -= frac[i];
		}
		x[frac_sort[N-1].i] = S;
	}
	
	return x;
}


/// The probability of the 	multinomal distribution
double multinomial_probability(unsigned int N_total, const vector <unsigned int> &x, const vector <double> &frac)
{
	auto sum = factorial(N_total);
	
	auto N = frac.size();
	for(auto i = 0u; i < N; i++){
		auto n = x[i];
		if(n > 0) sum += n*log(frac[i]);
		sum -= factorial(n);
	}
	
	return sum;
}

/// The probability of the 	multinomal distribution (without factorials)
double multinomial_samp_probability(const vector <unsigned int> &x, const vector <double> &frac)
{
	auto sum = 0.0;
	
	auto N = frac.size();
	for(auto i = 0u; i < N; i++){
		auto n = x[i];
		if(n > 0) sum += n*log(frac[i]);
	}
	
	return sum;
}


/// Works out which event occurs before the given one for a classification
unsigned int get_event_before(unsigned int e, const vector <Event> &ev)
{
	if(e == 0) return UNSET;
	auto cl = ev[e].cl;
	if(cl == UNSET) emsg("CL prob");
	e--;
	
	while(e > 0){
		auto &eve = ev[e];
		if(eve.cl == cl && (eve.type == M_TRANS_EV || eve.type == NM_TRANS_EV)) break;
		e--;
	}
	
	return e;
}


/// Works out which event occurs after the given one for a classification
unsigned int get_event_after(unsigned int e, const vector <Event> &ev)
{
	auto cl = ev[e].cl;

	e++;

	while(e < ev.size() && !((ev[e].type == M_TRANS_EV || ev[e].type == NM_TRANS_EV)
		&& ev[e].cl == cl) && ev[e].c_after != UNSET){
		e++;
	}
	if(e == ev.size()) return UNSET;
	return e;
}


/// Gets the log of a number but truncates if out of range
double get_log_zero_one(double val)
{
	if(val <= 0) return LI_WRONG;
	if(val >= 1) return 0;
	return log(val);
}


/// Samples from a list of posibilities
Poss sample_possibility(const vector <Poss> &pos)
{
	auto N = pos.size();
	
	if(N == 0) emsg("Cannot sample from possibility");
	
	vector <double> sum_store(N);
	
	auto sum = 0.0;
	for(auto j = 0u; j < N; j++){
		sum += pos[j].weight;
		sum_store[j] = sum;
	}
	if(sum == 0) emsg("Cannot sample from possibility");
	
	auto z = ran()*sum;
	
	auto k = 0u; while(k < N && z > sum_store[k]) k++;
	if(k == N) emsg("Cannot sample from possibility");	
	
	return pos[k];
}

/// Gets the log of the probability of observing 
double sample_probability(unsigned int i, const vector <Poss> &pos)
{
	auto sum = 0.0;
	auto jst = UNSET;
	for(auto j = 0u; j < pos.size(); j++){
		sum += pos[j].weight;
		if(pos[j].i == i) jst = j;
	}
	if(jst == UNSET) emsg("Could not find samp");
	
	return log(pos[jst].weight/sum);
}


/// Trncates a string to 40 characters
string trunc(const string &str)
{
	if(str.length() > TRUNC_MAX) return str.substr(0,TRUNC_MAX)+"...";
	return str;
}


/// CHecks that the branch probability is in range
void check_bp(double &bp)
{
	if(bp < 0 || bp > 1) emsg("Branching probability outside of range");
}						


/// Gets the likelihood for a transition
double nm_trans_like(TransType type, double dt, const vector <double> &ref_val)  
{		
	switch(type){
	case EXP_RATE: emsg("Should not be in NM"); return 0;
	case EXP_RATE_NM:	return exp_rate_probability(dt,ref_val[0]); 
	case GAMMA: return gamma_probability(dt,ref_val[0],ref_val[1]);
	case ERLANG: return gamma_probability(dt,ref_val[0],sqrt(1.0/ref_val[1])); 
	case LOG_NORMAL: return lognormal_probability(dt,ref_val[0],ref_val[1]);
	case WEIBULL:	return weibull_probability(dt,ref_val[0],ref_val[1]);
	case PERIOD: return period_probability(dt,ref_val[0]); 
	}
	return UNSET;
}

/// The log probability of sampling from a distribution
double prior_probability(double x, const Prior &pri, const vector <double> &param_val,const vector <Equation> &eqn)
{				
	switch(pri.type){
	case UNIFORM_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x < min || x > max) return -LARGE;
			return log(1.0/(max-min));
		}
		break;
		
	case EXP_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x < 0 || mean <= 0) return -LARGE;
			return exp_mean_probability(x,mean);
		}
		break;
	
	case NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto sd = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(sd <= 0) return -LARGE;
			return normal_probability(x,mean,sd);
		}
		break;
		
	case GAMMA_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || mean <= 0 || cv <= 0) return -LARGE;
			return gamma_probability(x,mean,cv);
		}
		break;
	
	case LOG_NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || mean <= 0 || cv <= 0) return -LARGE;
			return lognormal_probability(x,mean,cv);
		}
		break;
		
	case BETA_PR:
		{	
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			auto beta = eqn[pri.dist_param[1].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || x >= 1 || alpha <= 0 || beta <= 0) return -LARGE;
			return beta_probability(x,alpha,beta);
		}
		break;
	
	case BERNOULLI_PR:
		{	
			auto z = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x != 0 && x != 1) return -LARGE;
			if(z < 0 || z > 1) return -LARGE;
			return bernoulli_probability(x,z);
		}
		break;
		
	case FIX_PR:
		{
			auto val = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x != val) return -LARGE;
			return 0;
		}
		break;
		
	case FLAT_PR:
		{
			if(x <= 0) return -LARGE;
			return gamma_alpha_probability(x,1);
		}
		break;
		
	case DIRICHLET_PR:
		{
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param_only(param_val);
			if(x <= 0 || alpha <= 0) return -LARGE; 
			return gamma_alpha_probability(x,alpha);
		}
		break;
	}
	
	return UNSET;
}


/// Comma splits text but allows for quotation marks
vector <string> comma_split(string te)
{
	vector <string> vec;
	vec.push_back("");
	
	auto quote = 0u;
	auto num = 0u;
	for(auto i = 0u; i < te.length(); i++){
		auto cha = te.substr(i,1);
		if(cha == ","){
			if(quote == 0){ num++; vec.push_back("");}
			else vec[num] += cha;
		}
		else{
			if(cha == "\""){ 	
				quote = 1-quote;
			}
			else{
				vec[num] += cha;
			}
		}
	}

	return vec;
}


/// Adds like li1 to li2
void add_on_like(const Like li1, Like &li2)
{
	li2.init_cond += li1.init_cond;
	li2.init_cond_prior += li1.init_cond_prior;
	li2.obs += li1.obs;
	li2.prior += li1.prior; 
	li2.spline_prior += li1.spline_prior; 
	li2.dist += li1.dist;
	li2.markov += li1.markov; 
	li2.nm_trans += li1.nm_trans;
	li2.genetic_process += li1.genetic_process;
	li2.genetic_obs += li1.genetic_obs;
	li2.ie += li1.ie;
}


/// Calculates acceptance probability base don changes in likelihoods
double calc_al(const Like &like_ch, double dprob, const BurnInfo &burn_info)
{						
	if(false){
		return exp(like_ch.init_cond 
						+ like_ch.init_cond_prior
						+ like_ch.obs 
						+ like_ch.prior 
						+ like_ch.spline_prior 
						+ like_ch.dist 
						+ like_ch.markov 
						+ like_ch.nm_trans 
						+ like_ch.ie
						+ like_ch.genetic_process 
						+ like_ch.genetic_obs 
						+ dprob);
	}
	
	if(burn_info.dprob_suppress == true) if(dprob < -10) dprob = -10;
	
	return exp(burn_info.init_cond*like_ch.init_cond 
					+ burn_info.init_cond_prior*like_ch.init_cond_prior
					+ burn_info.obs*like_ch.obs 
					+ burn_info.prior*like_ch.prior 
					+ burn_info.spline_prior*like_ch.spline_prior 
					+ burn_info.dist*like_ch.dist 
					+ burn_info.markov*like_ch.markov 
					+ burn_info.nm_trans*like_ch.nm_trans 
					+ burn_info.ie*like_ch.ie
					+ burn_info.genetic_process*like_ch.genetic_process 
					+ burn_info.genetic_obs*like_ch.genetic_obs 
					+ dprob); 
}

/// Prints likelihood to the terminal
void print_like(const Like &like)
{
	cout << "Likelihood: ";
	if(like.init_cond != 0) cout << "ind_cond " << like.init_cond << "  ";
	if(like.init_cond_prior != 0) cout << "ind_cond_prior " << like.init_cond_prior << "  ";
	if(like.obs != 0) cout << "obs " << like.obs << "  ";
	if(like.prior != 0) cout << "prior " << like.prior << "  ";
	if(like.spline_prior != 0) cout << "spline_prior " << like.spline_prior << "  ";
	if(like.dist != 0) cout << "dist " << like.dist << "  ";
	if(like.markov != 0) cout << "markov " << like.markov << "  ";
	if(like.nm_trans != 0) cout << "nm_trans " << like.nm_trans << "  ";
	if(like.genetic_process  != 0) cout << "genetic_process  " << like.genetic_process  << "  ";
	if(like.genetic_obs  != 0) cout << "genetic_obs  " << like.genetic_obs  << "  ";
	if(like.ie  != 0) cout << "ie  " << like.ie  << "  ";
	cout << endl;
}


/// Prints a statement to terminal
void print(string te)
{
	if(!com_op && op()) cout << te << endl;
}


/// Prints a diagnostic statement to terminal
void print_diag(string te)
{
	if(print_diag_on && !com_op && op()) cout << te << endl;
}


/// Prints a diagnostic statement to terminal
void print(double num, string te)
{
	if(!com_op && op()) cout << num << " " << te << endl;
}


/// Gets the step use for printing sample
unsigned int get_op_step(unsigned int number)
{
	auto target = number/50;
	auto num = 1u;
	do{
		if(num > target) return num;
		if(2*num > target) return 2*num;
		if(5*num > target) return 5*num;
		num *= 10;
	}while(true);
}


/// Makes a list of a value which are true
vector <unsigned int> true_list(const vector <bool> &vec)
{
	vector <unsigned int> list;
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == true) list.push_back(i);
	}
	
	return list;
}


/// Makes a list of a value which are false
vector <unsigned int> false_list(const vector <bool> &vec)
{
	vector <unsigned int> list;
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i] == false) list.push_back(i);
	}
	
	return list;
}


/// Combines together a list based on two other (making sure not the same)
vector <unsigned int> combine(const vector <unsigned int> &vec1, const vector <unsigned int> &vec2)
{
	auto max = 0u;
	for(auto va : vec1) if(va > max) max = va;
	max++;
	vector <bool> map(max,false);
	for(auto va : vec1) map[va] = true;
	
	auto vec = vec1; 
	for(auto va : vec2){
		if(va >= max || map[va] == false) vec.push_back(va);
	}
	
	return vec;
}


/// Determines if a string begins with another string
bool begin_str(string st, string st2)
{
	auto tr = trim(st);
	if(st.length() < st2.length()) return false;
	if(st.substr(0,st2.length()) == st2) return true;
	return false;
}


/// Determines if a string ends with another string
bool end_str(string st, string st2)
{
	auto tr = trim(st);
	if(st.length() < st2.length()) return false;
	if(st.substr(st.length()-st2.length(),st2.length()) == st2) return true;
	return false;
}
