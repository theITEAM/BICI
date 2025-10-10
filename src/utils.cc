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
#ifdef MAC
#include <sys/time.h>
#include <sys/resource.h>
#endif

using namespace std;

#include "const.hh"
#include "utils.hh"


#ifdef WINDOWS
#include "windows.h"
#include "psapi.h"
#endif

#ifdef LINUX
#include "sys/types.h"
#include "sys/sysinfo.h"
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

default_random_engine generator;
std::uniform_real_distribution<> dist(0,1);
std::random_device rd; 
std::mt19937 gen(rd());

vector <double> log_sum;
vector <double> log_integer;

double percent_done;

//ofstream progress;
//bool progress_on = false;

/// Displays an error message which comes from an internal problem and kills otehr cores
void emsg(const string &msg)
{
	if(false) throw(std::runtime_error(msg));
	
	cout << endl;
	
	display_error(msg,true);
	if(false) raise(SIGABRT);
	
#ifdef USE_MPI
	//MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	int num;
	MPI_Comm_rank(MPI_COMM_WORLD,&num);
	cout << "Error from core: " << num << endl;	
	//MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
	exit (EXIT_FAILURE);
	//MPI_Finalize();
#else
	exit (EXIT_FAILURE);
#endif
}


/// Displays a human readable error message (this only exits when all core reach it)
void alert_input(const string &msg)
{
	if(false) throw(std::runtime_error(msg));
	if(op()) display_error(msg,false);
	if(false) raise(SIGABRT);
	
	end_code();
}


/// Displays a human readable error message and kill s other mpi process
void run_error(const string &msg)
{
	if(false) throw(std::runtime_error(msg));
	if(op()) display_error(msg,false);
	if(false) raise(SIGABRT);

	exit (EXIT_FAILURE);
}


/// Displays an error message for internal problem (this only exits when all core reach it)
void emsg_input(const string &msg)
{
	if(false) throw(std::runtime_error(msg));
	if(op()) display_error(msg,true);
	if(false) raise(SIGABRT);
	
	end_code();
}


/// Ends the code
void end_code()
{
	#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);  
  MPI_Finalize();
	exit (EXIT_SUCCESS);
#else
	exit (EXIT_SUCCESS);
#endif
}


/// Displays an error message
void display_error(string msg, bool internal)
{
	add_full_stop(msg);
	if(internal){
		msg = " This is an internal error to BICI and not something you have done wrong! Please send the BICI script to Chris for diagnosis. Error generated: "+msg;
	}
	
	if(com_op == true){
		cout << "<<ERROR>>" << endl <<  msg << endl << "<<END>>" << endl;
	}
	else{
		cout << endl << "\033[31m" << "ERROR:" << endl;
		cout << "\033[0m";
		cout << add_escape_char(msg);
		cout << endl;
	}
}


/// Display a warning message
void display_warning(string msg)
{
	add_full_stop(msg);
	
	if(com_op == true){
	}
	else{
		cout << "\033[35m";
		cout << "WARNING: ";
		cout << "\033[0m";
		cout << add_escape_char(msg);
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
vector<string> split(const string &s, char delimiter, bool notrim)
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
	
	if(notrim != true){
		for(auto &spl : splits) spl = trim(spl);
	}
	
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


/// Finds the index of a value in a vector
unsigned int find_in(const vector <Warn> &vec, const unsigned int val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i].th == val) return i;
	}
	return UNSET;
}


/// Finds the index of a value in a vector
unsigned int find_in(const vector <IEGref> &vec, const IEGref val)
{
	for(auto i = 0u; i < vec.size(); i++){
		if(vec[i].p == val.p && vec[i].i == val.i) return i;
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
	if(val == 0 && trim(st) != "0" && trim(st).substr(0,3) != "0.0" && trim(st).substr(0,3) != "0.") return UNSET;
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
void set_seed(const int core, const Details &details, unsigned int seed_tag)
{
	auto seed = details.seed;
	if(seed == UNSET) seed = SEED_DEFAULT;

	// Seed from tag loading BICI overides seed from file
	if(seed_tag != UNSET) seed = seed_tag;
	
	seed += 10000*core;

	generator.seed(seed);
	srand(seed);
}


/// Draws a random number between 0 and 1
double ran()
{
	double val;
	do{
		val = double(rand())/RAND_MAX;
		//val = dist(gen);
	}while(val == 0 || val == 1);

	return val;
}


/// Draws a normally distributed number with mean mu and standard deviation sd
double normal_sample(const double mean, const double sd, string &warn) 
{
	if(check_thresh(NORM_TE,SD_QU,sd,warn)) return UNSET;
	if(check_thresh(NORM_TE,NORM_MEAN_QU,mean,warn)) return UNSET;

	normal_distribution<double> distribution(mean, sd);
	return distribution(generator);
}


/// Returns an integer samping distribution with a normal envelope
int normal_int_sample(double si, string &warn)
{
	auto val = normal_sample(0,si,warn);
	if(val == UNSET) return UNSET;
	return int(val+0.5+LARGE_INT)-LARGE_INT;
}


/// The log of the probability from the normal distribution
double normal_probability(const double x, const double mean, const double sd)
{
	if(sd < SD_MIN) return LI_WRONG;
	auto var = sd*sd;
  return -0.5*log(2*MM_PI*var) - (x-mean)*(x-mean)/(2*var);
}


/*
/// The log of the probability from the normal distribution
double mvn_probability(const vector <double> x, const vector <double> mean, const vector < vector <double> > &covar)
{
	auto determinant_fast(const vector < vector <double> > &a)
	auto var = sd*sd;
	auto d = x.size();
	auto sum = 0.0;
	for(auto i = 0u; i < d; i++){
		for(auto j = 0u; j < d; j++){
			sum += (x[j]-mean[j])*covar[j][i]*(x[i]-mean[i]);
		}
	}
  return -(d/2)*log(2*MM_PI) - (d/2)*determinant_fast(covar) - 0.5*sum;
}
*/


/// Draws a log-normally distributed number with mean mu and standard deviation sd
double lognormal_sample(const double mean, const double cv, string &warn) 
{
	if(check_thresh(LOGNORM_TE,CV_QU,cv,warn)) return UNSET;
	if(check_thresh(LOGNORM_TE,MEAN_QU,mean,warn)) return UNSET;

	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	normal_distribution<double> distribution(mu, sqrt(var));
	auto val = exp(distribution(generator));
	if(val < TINY) return TINY;
	return val;
}


/// Probability of log-normally distribution
double lognormal_probability(const double x, const double mean, const double cv) 
{
	if(mean < MEAN_MIN || cv < CV_MIN || cv > CV_MAX || x < TINY) return LI_WRONG;
	
	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	auto logx = log(x);
	return  -0.5*log(2*MM_PI*var*x*x) - (logx-mu)*(logx-mu)/(2*var);
}


/// The lognormal probability of being xmin or above
double lognormal_upper_probability(const double xmin, const double mean, const double cv)
{
	if(mean < MEAN_MIN || cv < CV_MIN || cv > CV_MAX || xmin < TINY) return LI_WRONG;
	
	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	auto val = 1-0.5*erfc(-(log(xmin)-mu)/(sqrt(2*var)));
	if(std::isnan(val) || val < TINY) return LI_WRONG;
	return log(val);
}


/// The lognormal probability of being xmin or above (with no log taken)
double lognormal_upper_probability_no_log(const double xmin, const double mean, const double cv)
{
	if(mean < MEAN_MIN || mean > MEAN_MAX || cv < CV_MIN || cv > CV_MAX) return TINY;
	
	auto var = log(1+cv*cv);                // Works out variables on log scale
	auto mu = log(mean)-var/2;
	
	auto val = 1-0.5*erfc(-(log(xmin)-mu)/(sqrt(2*var)));
	if(std::isnan(val) || val < TINY) return TINY;
	return val;
}


/// Draws a Weibull distributed number with shape and scale parameters
double weibull_sample(const double scale, const double shape, string &warn) 
{
	if(check_thresh(WEIBULL_TE,SCALE_QU,scale,warn)) return UNSET;
	if(check_thresh(WEIBULL_TE,SHAPE_QU,shape,warn)) return UNSET;

	weibull_distribution<double> distribution(shape, scale);
	auto val = distribution(generator);
	if(val < TINY) return TINY;
	return val;
}


/// Probability of weibull sample
double weibull_probability(const double x, const double scale, const double shape) 
{
	if(shape < SHAPE_MIN || shape > SHAPE_MAX || 
	scale < SCALE_MIN || scale > SCALE_MAX) return LI_WRONG;
	
	return log(shape/scale)+(shape-1)*log(x/scale) - pow(x/scale,shape);
}


/// Probability of weibull sample beign xmin or above
double weibull_upper_probability(const double xmin, const double scale, const double shape) 
{
	if(scale < SCALE_MIN || scale > SCALE_MAX || shape < SHAPE_MIN || shape > SHAPE_MAX) return  LI_WRONG;
	
	auto val = -pow((xmin/scale),shape);
	if(std::isnan(val) || val < LI_WRONG) return LI_WRONG;
	return val;
}


/// Probability of weibull sample beign xmin or above (with no log taken)
double weibull_upper_probability_no_log(const double xmin, const double scale, const double shape) 
{
	if(scale < SCALE_MIN || scale > SCALE_MAX || shape < SHAPE_MIN || shape > SHAPE_MAX) return TINY;
	
	auto val = exp(-pow((xmin/scale),shape));
	if(std::isnan(val) || val < TINY) return TINY;
	return val;
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-b*x)
double gamma_sample(const double mean, const double cv, string &warn)
{
	if(check_thresh(GAMMA_TE,CV_QU,cv,warn)) return UNSET;
	if(check_thresh(GAMMA_TE,MEAN_QU,mean,warn)) return UNSET;

	gamma_distribution<double> distribution(1.0/(cv*cv),mean*cv*cv);
	auto val = distribution(generator);
	if(val < TINY) return TINY;
	return val;
}


/// The log of the probability from the gamma distribution
double gamma_probability(const double x, const double mean, const double cv)
{
	if(x < MEAN_MIN || mean < MEAN_MIN || cv < CV_MIN || cv > CV_MAX){
		return LI_WRONG;
	}
	
	auto shape = 1.0/(cv*cv);
	auto b = shape/mean;
	return (shape-1)*log(x) - b*x + shape*log(b) - lgamma(shape);
}


/// The gamma probability of being x or above
double gamma_upper_probability(const double xmin, const double mean, const double cv)
{
	if(xmin < MEAN_MIN || mean < MEAN_MIN || cv < CV_MIN || cv > CV_MAX) return LI_WRONG;
	
	auto shape = 1.0/(cv*cv);
	auto theta = mean/shape;
	auto value = boost::math::gamma_q(shape,xmin/theta);	
	if(std::isnan(value) || value < TINY) return LI_WRONG;
	return log(value);
}


/// The gamma probability of being x or above (with no log taken)
double gamma_upper_probability_no_log(const double xmin, const double mean, const double cv)
{
	if(xmin < MEAN_MIN || mean < MEAN_MIN || cv < CV_MIN || cv > CV_MAX) return TINY;
	
	auto shape = 1.0/(cv*cv);
	auto theta = mean/shape;
	auto val = boost::math::gamma_q(shape,xmin/theta);
	if(std::isnan(val) || val < TINY) return TINY;
	return val;
}


/// Draws a sample from the beta distribution x^(a-1)*(1-x)^(b-1)
double beta_sample(const double alpha, const double beta, string &warn)
{
	if(check_thresh(BETA_TE,ALPHA_QU,alpha,warn)) return UNSET;
	if(check_thresh(BETA_TE,BETA_QU,beta,warn)) return UNSET;

	gamma_distribution<double> distributionX(alpha,1), distributionY(beta,1);
	auto x = distributionX(generator);
	auto y = distributionY(generator);

	return x/(x+y);
}


/// The log of the probability from the beta distribution
double beta_probability(const double x, const double alpha, const double beta)
{
	if(x < TINY || x > 1-TINY) return LI_WRONG;
	if(alpha < ALPBETA_MIN || beta < ALPBETA_MIN ||
	   alpha > ALPBETA_MAX || beta > ALPBETA_MAX) return LI_WRONG;
	
  return (alpha-1)*log(x) + (beta-1)*log(1-x) + lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta);
}


/// Generates a sample from the Poisson distribution
unsigned int poisson_sample(const double lam, string &warn)
{
	if(check_thresh(POIS_TE,POIS_QU,lam,warn)) return UNSET;
	
  poisson_distribution<int> distribution(lam);
	return distribution(generator);
}


/// The log probability of the poisson distribution
double poisson_probability(const int i, const double lam)
{
	if(lam < LAM_MIN || lam > LAM_MAX) return LI_WRONG;
	
  if(lam == 0){
    if(i == 0) return 0;
    else return LI_WRONG;
  }
  else{
    if(i < 0) return LI_WRONG;
    else{
			return i*log(lam)-lam - lgamma(i+1);
		}
  }
}


/// The probability for above imin
double poisson_upper_probability_no_log(const int imin, const double lam)
{
	return 1-boost::math::gamma_q(imin+1,lam);
}


/// The log probability of the negative binomial distribution
double neg_binomial_probability(const int k, const double mean, const double p)
{	
	if(mean < MEAN_MIN || mean > MEAN_MAX || p < P_MIN || p > P_MAX) return LI_WRONG;
	
	auto r = round_int(mean*p/(1-p));
	if(r < 1) r = 1;
	return choose(k+r-1,k) + k*log(1-p) + r*log(p);
}


/// Generates a sample from the Bernoulli distribution
unsigned int bernoulli_sample(const double p, string &warn)
{
	if(check_thresh(BERN_TE,BERNP_QU,p,warn)) return UNSET;
	
	if(ran() < p) return 1;
	return 0;
}


/// The probability of a Bernoulli sample
double bernoulli_probability(unsigned int x, const double p)
{
	if(x == 0){
		if(p > 1-TINY) return LI_WRONG;
		return log(1-p);
	}
	
	if(p < TINY) return LI_WRONG;
	return log(p);
}


/// Samples from the exponential distribution with specified rate
double exp_rate_sample(const double rate, string &warn)
{
	if(check_thresh(EXP_RATE_TE,RATE_QU,rate,warn)) return UNSET;
	
	auto val = -log(ran())/rate;
	if(val < TINY) return TINY;
	return val;
}


/// Probability of exponential distribution
double exp_rate_probability(const double x, const double rate)
{
	if(rate < RATE_MIN) return LI_WRONG;
	 
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
	auto val = exp(-rate*xmin);
	if(std::isnan(val) || val < TINY) return TINY;
	return val;
}


/// Samples from the exponential distribution with specified mean
double exp_mean_sample(const double mean, string &warn)
{
	if(check_thresh(EXP_MEAN_TE,EXP_MEAN_QU,mean,warn)) return UNSET;
	
	return -log(ran())*mean;
}


/// Probability of exponential distribution
double exp_mean_probability(const double x, const double mean)
{
	if(x < 0 || mean < EXP_MEAN_MIN) return LI_WRONG;
	return -log(mean) - x/mean;
}


/// Probability of exponential distribution above xmin
double exp_mean_upper_probability(const double xmin, const double mean)
{
	if(xmin < 0 || mean < EXP_MEAN_MIN) return LI_WRONG;
	
	return -xmin/mean;
}


/// Probability of exponential distribution above xmin (with no log taken)
double exp_mean_upper_probability_no_log(const double xmin, const double mean)
{
	if(xmin < 0 || mean < EXP_MEAN_MIN) return TINY;
	
	auto val = exp(-xmin/mean);
	if(std::isnan(val) || val < TINY) return TINY;
	return val;
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-x)
double gamma_alpha_sample(const double alpha, string &warn)
{
	if(check_thresh(GAMMA_TE,ALPHA_QU,alpha,warn)) return UNSET;
	
	//if(alpha < TINY) emsg("Alpha must be positive");
	gamma_distribution<double> distribution(alpha,1);
	return distribution(generator);
}


/// Draws a sample from the gamma distribution x^(a-1)*exp(-x)
double gamma_alpha_probability(const double x, const double alpha)
{
  if(x <= 0 || alpha < ALPBETA_MIN || alpha > ALPBETA_MAX) return LI_WRONG;
  return (alpha-1)*log(x) - x - lgamma(alpha);
}


/// A sample from the binomial distribution
unsigned int binomial_sample(double p, unsigned int n, string &warn)
{
	if(n == 0) return 0;
	
	if(n > LARGE){ warn = "Binomual cannot sample from negative"; return UNSET;} 
	if(p < -TINY){ warn = "For the binomial distribution the probability must be postive"; return UNSET;}
	if(p > 1){
		warn = "For the binomial distribution the probability must be one or less";
		return UNSET;
	}
	
	if(p < TINY) return 0;
	if(p > ALMOST_ONE) return n;
	
	binomial_distribution<int> distribution(n,p);
	return distribution(generator);
}


/// A sample from the binomial distribution
double binomial_probability(unsigned int num, double p, unsigned int n)
{	
	if(p < 0 || p > 1) return LI_WRONG;

	if(p == 0){
		if(num == 0) return 0;
		else return LI_WRONG;
	}
	else{
		if(p == 1){
			if(num == n) return 0;
			else return LI_WRONG;
		}
		else{
			if(num > n) emsg("Problem");

			return lgamma(n+1) - lgamma(num+1) - lgamma(n-num+1) + num*log(p) + (n-num)*log(1-p);
		}
	}
}


/// Probability of period sample
double period_sample(const double time, string &warn) 
{
	if(check_thresh(PERIOD_TE,TIME_QU,time,warn)) return UNSET;
	return time;
}


/// Probability of period sample (Uses a very sharp spike prior)
double period_probability(const double x, const double time) 
{
	if(time < TIME_MIN) return LI_WRONG;
	auto dd = x - time;
	if(dd > 0) return -dd*LARGE;
	else return dd*LARGE;
}


/// Probability of period sample (Uses a very sharp spike prior)
double period_upper_probability(const double x, const double time) 
{
	if(x < time) return 0; 
	return LI_WRONG;
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


/// Compares two boolean vector to see if they are equal
bool equal_vec(const vector <bool> &vec1, const vector <bool> &vec2)
{
	if(vec1.size() != vec2.size()) return false;
	
	for(auto i = 0u; i < vec1.size(); i++){
		if(vec1[i] != vec2[i]) return false;
	}
	
	return true;
}


/// Gets a list up to a number
vector <unsigned int> get_list(unsigned int num)
{
	vector <unsigned int> list;
	for(auto i = 0u; i < num; i++) list.push_back(i);
	
	return list;
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
	DepInfo di; di.warn = "";

	string dep;

	auto ist = i;
	if(te.substr(i,1) == "("){	
		while(i < te.length() && te.substr(i,1) != ")") i++;
		if(i == te.length()) di.warn = "Brackets do not match";

		dep = te.substr(ist+1,i-(ist+1));
	}
	else{
		while(i < te.length() && !includes(not_allowed,te.substr(i,1)) && !str_eq(te,i,sigma) && te.substr(i,1) != "[") i++;
		dep = te.substr(ist,i-ist);
	}

	if(dep != ""){
		auto spl = split(dep,',',true);

		for(auto &de : spl){
			auto sh = 0u; while(sh < de.length() && de.substr(sh,1) == " ") sh++; 
			di.ipos.push_back(ist+sh);
			ist += de.length() + 1;
	
			de = trim(de);
			if(de == "") di.warn = "Dependency is not set"; 
		}
		
		for(auto d = 0u; d < spl.size(); d++){
			for(auto dd = d+1; dd < spl.size(); dd++){
				if(spl[d] == spl[dd]) di.warn = "Cannot have repeated index '"+spl[d]+"'"; 
			}
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
				auto di = get_dep_info(te,i,"$}[<(;");
			
				if(di.warn != ""){
					res.warn = "An error occured: ";
					if(type == "pop"){
						res.warn = "Population '{"+te.substr(ist,di.iend-ist)+"}' has an error: ";
					}
					if(type == "param"){
						res.warn = "Parameter '"+te.substr(istore+1,i-istore-2)+"' has a misspecified dependencyL ";
					}
					res.warn += di.warn;
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


/// Determines 
bool in_text(const string &te, unsigned int i, string te2)
{
	if(i+te2.length() > te.length()) return false;
	if(te.substr(i,te2.length()) == te2) return true;
	return false;
}


/// Constructs a swap template
SwapResult swap_template(string te, const vector <DepConv> &dep_conv)
{
	SwapResult res; res.warn = ""; res.done.resize(dep_conv.size(),false);
	
	auto icopy = 0u;
	
	auto i = 0u;
	while(i < te.length()){
		while(i < te.length() && te.substr(i,1) != "%" && te.substr(i,1) != "{" && !in_text(te,i,"Σ")) i++;
		if(i < te.length()){
			auto istore = i; 

			string type = "";
			if(in_text(te,i,"Σ")){  // Changes in a sum max
				while(i < te.length() && te.substr(i,1) != "[" && te.substr(i,1) != "(") i++;
				if(i < te.length() && te.substr(i,1) == "["){
					type = "sum";
					i++;
				}
			}
			else{
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
			}
			
			if(type != ""){
				auto end = "$}[<(;"; if(type == "sum") end = ",";
				
				auto ist = i;
				auto di = get_dep_info(te,i,end);

				if(di.warn != ""){
					res.warn = "An error occured: ";
					if(type == "pop"){
						res.warn = "Population '{"+te.substr(ist,di.iend-ist)+"}' has an error: ";
					}
					if(type == "param"){
						res.warn = "Parameter '"+te.substr(istore+1,i-istore-2)+"' has a misspecified dependency: ";
					}
					if(type == "sum"){
						res.warn = "Sum '"+te.substr(istore+1,i-istore-2)+"' has a misspecified dependency: ";
					}
					res.warn += di.warn;
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
		//emsg("Split");
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


/// Adds a value to a vector (if it doesn't already exist)
unsigned int add_to_vec(vector <PopTransRef> &vec, unsigned int p, unsigned int tr)
{
	auto i = 0u; 
	while(i < vec.size() && !(vec[i].p == p && vec[i].tr == tr)) i++;
	if(i == vec.size()){
		PopTransRef pref; pref.p = p; pref.tr = tr; 
		vec.push_back(pref);
	}	
	return i;
}


/// Adds a value to a vector (if it doesn't already exist)
unsigned int add_to_vec(vector <PopMarkovEqnRef> &vec, unsigned int p, unsigned int e)
{
	auto i = 0u; 
	while(i < vec.size() && !(vec[i].p == p && vec[i].e == e)) i++;
	if(i == vec.size()){
		PopMarkovEqnRef pref; pref.p = p; pref.e = e; 
		vec.push_back(pref);
	}	
	return i;
}


/// Adds a value to a vector (if it doesn't already exist)
unsigned int add_to_vec(vector <IEGref> &vec, const IEGref &val)
{
	auto i = 0u; 
	while(i < vec.size() && !(vec[i].p == val.p && vec[i].i == val.i)) i++;
	if(i == vec.size()){
		vec.push_back(val);
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
bool dif2(double a, double b, double thresh)
{
	if(std::isnan(a)){ return true;}
	if(std::isnan(b)){ return true;}
	
	auto dif = a-b; if(dif < 0) dif = -dif;
	
	if(dif > thresh){
		if(a < 0) a = -a;
		if(b < 0) b = -b;
		
		auto frac = dif/(a+b);
		if(frac > 0.00000001) return true;
	}
	return false;
}


/// Calculates if two numbers is different (subject to numerical noise) 
bool dif(double a, double b, double thresh)
{
	if(std::isnan(a)) return true;
	if(std::isnan(b)) return true;
	
	auto dif = a-b; if(dif < 0) dif = -dif;
	
	if(dif > thresh){
		if(a < 0) a = -a;
		if(b < 0) b = -b;
		
		auto frac = dif/(a+b);
		if(dif == DIF_THRESH){
			if(frac > 0.00000001) return true;
		}
		else{
			if(frac > 0.00001) return true;
		}
	}
	return false;
}


/// Calculates if two numbers is different (subject to numerical noise) 
bool dif(const vector <double> &a, const vector <double> &b, double thresh)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(dif(a[j],b[j],thresh)) return true;
	}
	return false;
}


/// returns which is different
unsigned int which_dif(const vector <double> &a, const vector <double> &b, double thresh)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(dif(a[j],b[j],thresh)) return j;
	}
	return UNSET;
}


/// Calculates if two numbers is different (subject to numerical noise) 
bool dif(const vector < vector <double> > &a, const vector < vector <double> > &b, double thresh)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(a[j].size() != b[j].size()) return true;
		for(auto i = 0u; i < a[j].size(); i++){
			if(dif(a[j][i],b[j][i],thresh)) return true;
		}
	}
	return false;
}


/// Calculates if two numbers is different (subject to numerical noise) 
vector <unsigned int> which_dif(const vector < vector <double> > &a, const vector < vector <double> > &b, double thresh)
{
	vector <unsigned int> vec;
	
	if(a.size() != b.size()) return vec;
	for(auto j = 0u; j < a.size(); j++){
		if(a[j].size() != b[j].size()) return vec;
		for(auto i = 0u; i < a[j].size(); i++){
			if(dif(a[j][i],b[j][i],thresh)){
				vec.push_back(j); vec.push_back(i); 
				return vec;
			}
		}
	}
	return vec;
}


/// Calculates if two unsigned ints are different 
bool difi(unsigned int a, unsigned int b)
{
	if(a != b) return true;
	return false;
}


/// Calculates if two vectors are different 
bool difi(const vector <unsigned int> &a, const vector <unsigned int> &b)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(difi(a[j],b[j])) return true;
	}
	return false;
}


/// Calculates if two numbers is different 
bool difi(const vector < vector <unsigned int> > &a, const vector < vector <unsigned int> > &b)
{
	if(a.size() != b.size()) return true;
	for(auto j = 0u; j < a.size(); j++){
		if(a[j].size() != b[j].size()) return true;
		for(auto i = 0u; i < a[j].size(); i++){
			if(difi(a[j][i],b[j][i])) return true;
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
		if(st.substr(i,st1.length()) == st1){
			st = st.substr(0,i)+st2+st.substr(i+st1.length());
			i += st2.length();
		}
		else i++;
	}
	return st;
}

/// Replaces any tab characters with spaces
void remove_tab(string &st)
{
	st = replace(st,"\t"," ");
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
	case MARKOV_POP_NOPOP_AFFECT: return "MARKOV_POP_NOPOP_AFFECT";
	case MARKOV_POP_LINEAR_AFFECT: return "MARKOV_POP_LINEAR_AFFECT";
	case LIKE_IE_AFFECT: return "LIKE_IE_AFFECT";
	case EXP_IE_AFFECT: return "EXP_IE_AFFECT";
	case OMEGA_AFFECT: return "OMEGA_AFFECT";
	case POP_AFFECT: return "POP_AFFECT";
	case SPLINE_PRIOR_AFFECT: return "SPLINE_PRIOR_AFFECT";
	case IEG_PRIOR_AFFECT: return "IEG PRIOR_AFFECT";
	case PRIOR_AFFECT: return "PRIOR_AFFECT";
	case DIST_AFFECT: return "DIST_AFFECT";
	case EXP_FE_AFFECT: return "EXP_FE_AFFECT";
	case INDFAC_INT_AFFECT: return "INDFAC_INT_AFFECT";
	case DIV_VALUE_AFFECT: return "DIV_VALUE_AFFECT";
	case DIV_VALUE_NOPOP_AFFECT: return "DIV_VALUE_NOPOP_AFFECT";
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
	case DIV_VALUE_AFFECT: case DIV_VALUE_NOPOP_AFFECT: 
	case DIV_VALUE_LINEAR_AFFECT: 
	case MARKOV_LIKE_AFFECT: case POP_AFFECT:
	case MARKOV_POP_AFFECT: case MARKOV_POP_NOPOP_AFFECT: case MARKOV_POP_LINEAR_AFFECT:
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


/// Finds the mean value for an array
double mean(const vector <double> &val)
{
	auto sum = 0.0;
	for(auto v : val) sum += v;
	return sum/val.size();
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


/// Rounds a number down to the nearest integer 
int floor_int(double val)
{
	return int(val+LARGE_INT)-LARGE_INT;
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
string tstr(double value, unsigned int dig)
{
	stringstream ss; 
	if(dig != UNSET) ss << setprecision(dig); 
	ss << value; 
	return ss.str();
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


/// Determines the core number
unsigned int core()
{
#ifdef USE_MPI
	int num;
	MPI_Comm_rank(MPI_COMM_WORLD,&num); 
	return (unsigned int) num;
#else		
	return 0;
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


/// Constructs a sequention vector, e.g. n1,2,3,...n2-1
vector <unsigned int> seq_vec(unsigned int n1, unsigned int n2)
{
	vector <unsigned int> vec;
	for(auto i = n1; i < n2; i++) vec.push_back(i);
	
	return vec;
}


/// Constructs a vector with all one
vector <bool> true_vec(unsigned int N)
{
	vector <bool> vec(N,true);
	return vec;
}


/// Converts a text string to a prior specification
Prior convert_text_to_prior(string te, unsigned int line_num, string in, bool dist)
{
	Prior pri;
	te = trim(te);
	pri.name = te;
	pri.in = in;
	
	auto spl = split(te,'(');
	if(spl.size() != 2){ pri.error = "Syntax error"; return pri;}
	
	if(spl[1].size() == 0){ pri.error = "Syntax error"; return pri;}
	if(spl[1].substr(spl[1].length()-1,1) != ")"){ pri.error = "Syntax error"; return pri;}
	
	auto type = toLower(spl[0]);
	if(type == "dir") type = "dirichlet";
	if(type == "bern") type = "bernoulli";
	
	auto fl = false;
	if(type == "uniform"){ pri.type = UNIFORM_PR; fl = true;}
	if(type == "inverse"){ pri.type = INVERSE_PR; fl = true;}
	if(type == "power"){ pri.type = POWER_PR; fl = true;}
	if(type == "mvn-jeffreys"){ pri.type = MVN_JEF_PR; fl = true;}
	if(type == "mvn-uniform"){ pri.type = MVN_UNIFORM_PR; fl = true;}
	if(type == "exp"){ pri.type = EXP_PR; fl = true;}
	if(type == "normal"){ pri.type = NORMAL_PR; fl = true;}
	if(type == "gamma"){ pri.type = GAMMA_PR; fl = true;}
	if(type == "log-normal"){ pri.type = LOG_NORMAL_PR; fl = true;}
	if(type == "beta"){ pri.type = BETA_PR; fl = true;}
	if(type == "bernoulli"){ pri.type = BERNOULLI_PR; fl = true;}
	if(type == "fix"){ pri.type = FIX_PR; fl = true;}
	if(type == "dirichlet"){ pri.type = DIRICHLET_PR; fl = true;}
	if(type == "mdir"){ pri.type = MDIR_PR; fl = true;}
	
	if(fl == false){
		pri.error = "Distribution '"+type+"' not recognised"; return pri;
	}
	
	auto bra = spl[1].substr(0,spl[1].length()-1);
	auto spl2 = split(bra,',');

	switch(pri.type){
	case MVN_COR_PR: emsg("Correlation not defined"); break;
	case INVERSE_PR: case MVN_JEF_PR: case MVN_UNIFORM_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
			
			if(!dist){
				if(number(spl2[0]) == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				if(number(spl2[1]) == UNSET){ 
					pri.error = "The value '"+spl2[1]+"' must be a number"; return pri;
				}
				if(number(spl2[0]) <= 0){  
					pri.error = "Minimum must be a positive number"; return pri;
				}
				if(number(spl2[1]) <= 0){  
					pri.error = "Maximum must be a positive number"; return pri;
				}
				if(number(spl2[0]) >= number(spl2[1])){  
					pri.error = "Minimum must be smaller than maximum";return pri;
				}
			}
		}
		break;
		
	case POWER_PR:
		{
			if(spl2.size() != 3){ pri.error = "Expected three values in the brackets"; return pri;}
			
			if(!dist){
				if(number(spl2[0]) == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				if(number(spl2[1]) == UNSET){ 
					pri.error = "The value '"+spl2[1]+"' must be a number"; return pri;
				}
				if(number(spl2[2]) == UNSET){ 
					pri.error = "The value '"+spl2[2]+"' must be a number"; return pri;
				}
				if(number(spl2[0]) <= 0){  
					pri.error = "Minimum must be a positive number"; return pri;
				}
				if(number(spl2[1]) <= 0){  
					pri.error = "Maximum must be a positive number"; return pri;
				}
				if(number(spl2[0]) >= number(spl2[1])){  
					pri.error = "Minimum must be smaller than maximum";return pri;
				}
			}
		}
		break;
		
	case FIX_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
			
			if(!dist){
				if(number(spl2[0]) == UNSET) pri.error = "The value '"+spl2[0]+"' must be a number";
			}
		}
		break;
	
	case UNIFORM_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
			
			if(!dist){
				if(number(spl2[0]) == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				if(number(spl2[1]) == UNSET){ 
					pri.error = "The value '"+spl2[1]+"' must be a number"; return pri;
				}
				if(number(spl2[0]) >= number(spl2[1])){  
					pri.error = "Minimum must be smaller than maximum";return pri;
				}
			}
		}
		break;
	
	case EXP_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
			
			auto mean = number(spl2[0]);
			
			if(!dist){
				if(mean == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				
				if(mean <= 0){
					pri.error = "Mean must be positive"; return pri;
				}
			}
			
			if(mean != UNSET && check_thresh(EXP_MEAN_TE,EXP_MEAN_QU,mean,pri.error)) return pri; 		
		}
		break;
		
	case NORMAL_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
			
			auto mean = number(spl2[0]);
			auto sd = number(spl2[1]);
	
			if(!dist){
				if(mean == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				if(sd == UNSET){
					pri.error = "The value '"+spl2[1]+"' must be a number"; return pri;
				}
			
				if(sd <= 0){
					pri.error = "Standard deviation must be positive"; return pri;
				}
			}
			
			if(sd != UNSET && check_thresh(NORM_TE,SD_QU,sd,pri.error)) return pri;		
			if(mean != UNSET && check_thresh(NORM_TE,NORM_MEAN_QU,mean,pri.error)) return pri; 		
		}
		break;
		
	case GAMMA_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
			
			auto mean = number(spl2[0]);
			auto cv = number(spl2[1]);
	
			if(!dist){
				if(mean == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				if(cv == UNSET){
					pri.error = "The value '"+spl2[1]+"' must be a number"; return pri;
				}
			
				if(mean <= 0){
					pri.error = "Mean must be positive"; return pri;
				}
				
				if(cv <= 0){
					pri.error = "Standard deviation must be positive"; return pri;
				}
			}
			
			if(cv != UNSET && check_thresh(GAMMA_TE,CV_QU,cv,pri.error)) return pri; 		
			if(mean != UNSET && check_thresh(GAMMA_TE,MEAN_QU,mean,pri.error)) return pri; 		
		}
		break;
		
	case LOG_NORMAL_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
			
			auto mean = number(spl2[0]);
			auto cv = number(spl2[1]);
	
			if(!dist){
				if(mean == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				if(cv == UNSET){
					pri.error = "The value '"+spl2[1]+"' must be a number"; return pri;
				}
			
				if(mean <= 0){
					pri.error = "Mean must be positive"; return pri;
				}
				
				if(cv <= 0){
					pri.error = "Standard deviation must be positive"; return pri;
				}
			}
			
			if(cv != UNSET && check_thresh(LOGNORM_TE,CV_QU,cv,pri.error)) return pri; 			
			if(mean != UNSET && check_thresh(LOGNORM_TE,MEAN_QU,mean,pri.error)) return pri; 		
		}
		break;
		
	case BETA_PR:
		{
			if(spl2.size() != 2){ pri.error = "Expected two values in the brackets"; return pri;}
			
			auto alpha = number(spl2[0]);
			auto beta = number(spl2[1]);
	
			if(!dist){
				if(alpha == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				if(beta == UNSET){
					pri.error = "The value '"+spl2[1]+"' must be a number"; return pri;
				}
			
				if(alpha <= 0){
					pri.error = "Alpha must be positive"; return pri;
				}
			
				if(beta <= 0){
					pri.error = "Beta must be positive"; return pri;
				}
			}
			
			if(beta != UNSET && check_thresh(BETA_TE,BETA_QU,beta,pri.error)) return pri; 		
			if(alpha != UNSET && check_thresh(BETA_TE,ALPHA_QU,alpha,pri.error)) return pri; 		
		}
		break;

	case BERNOULLI_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
			
			auto mean = number(spl2[0]);
		
			if(!dist){
				if(mean == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
			
				if(mean < 0 || mean > 1){
					pri.error = "Mean must be between zero and one"; return pri;
				}
			}
			
			if(mean != UNSET && check_thresh(BERN_TE,BERNP_QU,mean,pri.error)) return pri; 		
		}
		break;
	
	case DIRICHLET_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
			
			auto alpha = number(spl2[0]);
		
			if(!dist){
				if(alpha == UNSET){
					pri.error = "The value '"+spl2[0]+"' must be a number"; return pri;
				}
				
				if(alpha <= 0){
					pri.error = "Alpha must be positive"; return pri;
				}
			}
		}
		break;
		
	case MDIR_PR:
		{
			if(spl2.size() != 1){ pri.error = "Expected one value in the brackets"; return pri;}
			
			if(number(spl2[0]) == UNSET) pri.error = "The value '"+spl2[0]+"' must be a number";
		}
		break;
		
	case UNSET_PR: emsg("Should not be unset"); break;
	}

	for(auto i = 0u; i < spl2.size(); i++){
		auto ei = add_equation_info(spl2[i],DIST,line_num);
		if(ei.error){
			pri.error = ei.emsg;
			return pri;
		}
		pri.dist_param.push_back(ei);
	}
	
	return pri;
}
	
	
/// Converts a prior to a string
string get_prior_string(Prior prior)
{
	switch(prior.type){
	case MVN_COR_PR: emsg("Prior not defined"); break;
	case MVN_JEF_PR:
		return "mvn-jeffreys("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";
	
	case MVN_UNIFORM_PR:
		return "mvn-uniform("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";
	
	case INVERSE_PR:
		return "inverse("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+")";

	case POWER_PR:
		return "power("+prior.dist_param[0].te_raw+","+prior.dist_param[1].te_raw+","+prior.dist_param[2].te_raw+")";

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
	
	case DIRICHLET_PR:
		return "dir("+prior.dist_param[0].te_raw+")";
		
	case MDIR_PR:
		return "mdir("+prior.dist_param[0].te_raw+")";	
		
	case UNSET_PR: 
		emsg("Should not be unset"); 
		return "";
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
	
	string warn;
	
	vector <double> vec(N);
	auto sum = 0.0;
	for(auto i = 0u; i < N; i++){
		auto al = alpha[i];
		if(al == ALPHA_ZERO) vec[i] = 0;
		else{
			auto val = gamma_alpha_sample(al,warn);
			if(val == UNSET) emsg("Dirichlet distribution sampling error: "+warn); 
			vec[i] = val;
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
bool Frac_ord (const FracSort &fs1, const FracSort &fs2)                      
{ return (fs1.frac > fs2.frac); };  


/// Samples from a multinomial probability distribution
vector <unsigned int>  multinomial_sample(unsigned int N_total, const vector <double> &frac)
{
	auto N = frac.size();
	
	vector <unsigned int> x(N,0);
	
	string warn;
	
	if(false){ // Direct way
		vector <double> frac_sum(N);
		auto sum = 0.0;
		for(auto i = 0u; i < N; i++){
			sum += frac[i];
			frac_sum[i] = sum;
		}
		if(dif(sum,1,DIF_THRESH)) emsg("sum not one");
		
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
			if(r > 0){
				n = binomial_sample(ALMOST_ONE*frac[i]/r,S,warn);
				if(n == UNSET) emsg("Problem with multinomial sampling: "+warn);
			}
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
double nm_trans_like(TransType type, double dtdiv, double dt, const vector <double> &ref_val)  
{		
	switch(type){
	case EXP_RATE: case EXP_MEAN: emsg("Should not be in NM"); return 0;
	case EXP_RATE_NM:	return exp_rate_probability(dtdiv,ref_val[0]*dt); 
	case EXP_MEAN_NM:	return exp_mean_probability(dtdiv,ref_val[0]/dt); 
	case GAMMA: return gamma_probability(dtdiv,ref_val[0]/dt,ref_val[1]);
	case ERLANG: return gamma_probability(dtdiv,ref_val[0]/dt,sqrt(1.0/ref_val[1])); 
	case LOG_NORMAL: return lognormal_probability(dtdiv,ref_val[0]/dt,ref_val[1]);
	case WEIBULL:	return weibull_probability(dtdiv,ref_val[0]/dt,ref_val[1]);
	case PERIOD: return period_probability(dtdiv,ref_val[0]/dt); 
	}
	return UNSET;
}

/// The log probability of sampling from a distribution
double prior_probability(double x, const Prior &pri, const vector <double> &precalc, const vector <Equation> &eqn)
{				
	switch(pri.type){
	case MVN_JEF_PR: case MVN_UNIFORM_PR: case MVN_COR_PR: // This is done collectively
		return 0;
		
	case INVERSE_PR: 
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(min <= 0 || max <= 0) return -LARGE;
			if(x < min || x > max) return -LARGE;
			return log(1.0/(x*log(max/min)));
		}
		break;
		
	case POWER_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			auto power = eqn[pri.dist_param[2].eq_ref].calculate_param(precalc);
			if(min <= 0 || max <= 0) return -LARGE;
			if(x < min || x > max) return -LARGE;
			if(power == -1) return log(1.0/(x*log(max/min)));
			else{
				auto c = (power+1)/(pow(max,power+1)-pow(min,power+1));
				return power*log(x)+log(c);
			}
		}
		break;
		
	case UNIFORM_PR:
		{	
			auto min = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto max = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(x < min || x > max) return -LARGE;
			return log(1.0/(max-min));
		}
		break;
		
	case EXP_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			if(x < 0 || mean <= 0) return -LARGE;
			return exp_mean_probability(x,mean);
		}
		break;
	
	case NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto sd = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(sd <= 0) return -LARGE;
			return normal_probability(x,mean,sd);
		}
		break;
		
	case GAMMA_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(x <= 0 || mean <= 0 || cv <= 0) return -LARGE;
			return gamma_probability(x,mean,cv);
		}
		break;
	
	case LOG_NORMAL_PR:
		{	
			auto mean = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto cv = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(x <= 0 || mean <= 0 || cv <= 0) return -LARGE;
			return lognormal_probability(x,mean,cv);
		}
		break;
		
	case BETA_PR:
		{	
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			auto beta = eqn[pri.dist_param[1].eq_ref].calculate_param(precalc);
			if(x <= 0 || x >= 1 || alpha <= 0 || beta <= 0) return -LARGE;
			return beta_probability(x,alpha,beta);
		}
		break;
	
	case BERNOULLI_PR:
		{	
			auto z = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			if(x != 0 && x != 1) return -LARGE;
			if(z < 0 || z > 1) return -LARGE;
			return bernoulli_probability(x,z);
		}
		break;
		
	case FIX_PR:
		{
			auto val = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			if(x != val) return -LARGE;
			return 0;
		}
		break;
		
	case DIRICHLET_PR:
		{
			auto alpha = eqn[pri.dist_param[0].eq_ref].calculate_param(precalc);
			if(x <= 0 || alpha <= 0) return -LARGE; 
			return gamma_alpha_probability(x,alpha);
		}
		break;
		
	case MDIR_PR:
		emsg("Should not be mdir");
		break;
		
	case UNSET_PR:
		emsg("Should not be unset"); 
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
	li2.prior_bounded += li1.prior_bounded; 
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
						+ like_ch.prior_bounded 
						+ like_ch.spline_prior 
						+ like_ch.dist 
						+ like_ch.markov 
						+ like_ch.nm_trans 
						+ like_ch.ie
						+ like_ch.genetic_process 
						+ like_ch.genetic_obs 
						+ dprob);
	}
	
	//if(burn_info.dprob_suppress == true) if(dprob < -10) dprob = -10;
	
	return exp(burn_info.init_cond*like_ch.init_cond 
					+ burn_info.init_cond_prior*like_ch.init_cond_prior
					+ burn_info.obs*like_ch.obs 
					+ burn_info.prior*like_ch.prior 
				  + burn_info.prior_bounded*like_ch.prior_bounded 
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
	if(like.prior_bounded != 0) cout << "prior_bounded " << like.prior_bounded << "  ";
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


/// Gets the core
unsigned int get_core()
{
#ifdef USE_MPI
	int num;
	MPI_Comm_rank(MPI_COMM_WORLD,&num);
	return num;
#endif
	return 0;
}


/// Prints a diagnostic statement to terminal
void print_diag(string te)
{
	if(print_diag_on && !com_op && op()) cout << te << endl;
	//cout << te << endl; 
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


/// Checks if value is within thresholds
bool check_thresh(DistText dist, DistQuant dq, double val, string &err)
{
	if(val == UNSET){ err = "Value is unset"; return true;}
	
	double min_th = UNSET, max_th = UNSET;
	
	switch(dq){
	case MEAN_QU: 
		if(val < MEAN_MIN) min_th = MEAN_MIN;
		if(val > MEAN_MAX) max_th = MEAN_MAX; 
		break;
		
	case NORM_MEAN_QU: 
		if(val < NORM_MEAN_MIN) min_th = NORM_MEAN_MIN;
		if(val > NORM_MEAN_MAX) max_th = NORM_MEAN_MAX; 
		break;
	
	case SD_QU: 
		if(val < SD_MIN) min_th = SD_MIN;
		if(val > SD_MAX) max_th = SD_MAX; 
		break;
	
	case CV_QU: 
		if(val < CV_MIN) min_th = CV_MIN;
		if(val > CV_MAX) max_th = CV_MAX; 
		break;
		
	case SHAPE_QU: 
		if(val < SHAPE_MIN) min_th = SHAPE_MIN;
		if(val > SHAPE_MAX) max_th = SHAPE_MAX; 
		break;
		
	case SCALE_QU: 
		if(val < SCALE_MIN) min_th = SCALE_MIN;
		if(val > SCALE_MAX) max_th = SCALE_MAX; 
		break;
	
	case ALPHA_QU: case BETA_QU:
		if(val < ALPBETA_MIN) min_th = ALPBETA_MIN;
		if(val > ALPBETA_MAX) max_th = ALPBETA_MAX; 
		break;
		
	case P_QU: 
		if(val < P_MIN) min_th = P_MIN;
		if(val > P_MAX) max_th = P_MAX; 
		break;
		
	case BERNP_QU: 
		if(val < 0) min_th = 0;
		if(val > 1) max_th = 1; 
		break;
	
	case RATE_QU: 
		if(val < RATE_MIN) min_th = RATE_MIN;
		break;
		
	case EXP_MEAN_QU: 
		if(val < EXP_MEAN_MIN) min_th = EXP_MEAN_MIN;
		break;
		
	case TIME_QU: 
		if(val < TIME_MIN) min_th = TIME_MIN;
		break;
		
	case POIS_QU: 
		if(val < LAM_MIN) min_th = LAM_MIN;
		if(val > LAM_MAX) max_th = LAM_MAX; 
		break;
	}

	string dist_te;
	if(min_th != UNSET || max_th != UNSET){
		switch(dist){
		case LOGNORM_TE: dist_te = "Log-normal"; break;
		case NORM_TE: dist_te = "Normal"; break;
		case WEIBULL_TE: dist_te = "Weibull"; break;
		case GAMMA_TE: dist_te = "Gamma"; break;
		case BETA_TE: dist_te = "Beta"; break;
		case NEGBINO_TE: dist_te = "Negitive-binomial"; break;
		case BERN_TE: dist_te = "Bernoulli"; break;
		case EXP_RATE_TE: dist_te = "Exponential"; break;
		case EXP_MEAN_TE: dist_te = "Exponential"; break;
		case POIS_TE: dist_te = "Poisson"; break;
		case PERIOD_TE: dist_te = "Period"; break;
		} 
		
		string quant;
		switch(dq){
		case MEAN_QU: quant = "mean"; break;
		case NORM_MEAN_QU: quant = "mean"; break;
		case SD_QU: quant = "standard deviation"; break;
		case CV_QU: quant = "coefficient of variation"; break;
		case SHAPE_QU: quant = "shape"; break;
		case SCALE_QU: quant = "scale"; break;
		case ALPHA_QU: quant = "alpha"; break;
		case BETA_QU: quant = "beta"; break;
		case P_QU: quant = "probability"; break;
		case BERNP_QU: quant = "probability"; break;
		case RATE_QU: quant = "rate"; break;
		case EXP_MEAN_QU: quant = "mean"; break;
		case POIS_QU: quant = "mean"; break;
		case TIME_QU: quant = "time"; break;
		}		
		
		if(min_th != UNSET){
			err = dist_te+" "+quant+" has the value '"+tstr(val)+"' which should not be below the threshold '"+tstr(min_th)+"'";
		}
		
		if(max_th != UNSET){
			err = dist_te+" "+quant+" has the value '"+tstr(val)+"' which should not be above the threshold '"+tstr(max_th)+"'";
		}
	}
	
	if(err != "") return true;
	return false;
}


/// Truncates a text string to a certain length
string trunc(string te, unsigned int len)
{
	if(te.length() < len) return te;
	return te.substr(0,len)+"...";
}


/// Returns a list of escape characters
vector< vector <string> > get_escape_char()
{
	vector< vector <string> > escape_char;
	
	escape_char.push_back({"\\alpha","α"});
	escape_char.push_back({"\\beta","β"});
	escape_char.push_back({"\\gamma","γ"});
	escape_char.push_back({"\\Gamma","Γ"});
	escape_char.push_back({"\\delta","δ"});
	escape_char.push_back({"\\Delta","Δ"});
	escape_char.push_back({"\\epsilon","ε"});
	escape_char.push_back({"\\zeta","ζ"});
	escape_char.push_back({"\\eta","η"});
	escape_char.push_back({"\\Eta","Η"});
	escape_char.push_back({"\\theta","θ"});
	escape_char.push_back({"\\Theta","Θ"});
	escape_char.push_back({"\\iota","ι"});
	escape_char.push_back({"\\kappa","κ"});
	escape_char.push_back({"\\lambda","λ"});
	escape_char.push_back({"\\Lambda","Λ"});
	escape_char.push_back({"\\mu","μ"});
	escape_char.push_back({"\\nu","ν"});
	escape_char.push_back({"\\xi","ξ"});
	escape_char.push_back({"\\Xi","Ξ"});
	escape_char.push_back({"\\omicron","ο"});
	escape_char.push_back({"\\pi","π"});
	escape_char.push_back({"\\Pi","Π"});
	escape_char.push_back({"\\rho","ρ"});
	escape_char.push_back({"\\sigma","σ"});
	//escape_char.push_back({"\\Sigma","Σ"});
	escape_char.push_back({"\\tau","τ"});
	escape_char.push_back({"\\upsilon","υ"});
	escape_char.push_back({"\\phi","φ"});
	escape_char.push_back({"\\Phi","Φ"});
	escape_char.push_back({"\\chi","χ"});
	escape_char.push_back({"\\psi","ψ"});
	escape_char.push_back({"\\Psi","Ψ"});
	escape_char.push_back({"\\omega","ω"});	
	escape_char.push_back({"\\Omega","Ω"});	
	escape_char.push_back({"\\sum","Σ"});		
	escape_char.push_back({"\\int","∫"});		
	
	return escape_char;
}


/// Replaces greek letters with escape characters
string add_escape_char(string te)
{
	auto escape_char = get_escape_char();

	for(auto i = 0u; i < escape_char.size(); i++){
		auto st1 = escape_char[i][0];
		auto st2 = escape_char[i][1];
		if(includes(te,st2)){
			te = replace(te,st2,st1);
		}
	
		te = replace(te,"×","*");
		te = replace(te,"→","->");
	}
	
	return te;
}


/// Replaces right arrow
string replace_arrow(string te)
{
	return replace(te,"→","->");
}


/*
/// Initialises a file for outputting progress
void start_progess_file(string file)
{
	progress_on = true;
	progress.open(file);	
}
*/



/// Resets the percentage
void percentage_start(PercentType type, unsigned int gen)
{
	if(!op()) return;
	
	if(com_op){
		switch(type){
		case LOAD_PER: cout << "<CREATING>" << endl; break;
		case INIT_PER: cout << "<INIT>" << endl; break;
		case RUN_PER: cout << "<RUNNING>"  << endl; break;
		case RUN_GEN_PER: cout << "<RUNGEN>" << gen << endl; break;
		case ANNEAL_PER: cout << "<ANNEALING>" << endl; break;
		case OUTPUT_PER: cout << "<OUTPUTTING>" << endl; break;
		}
	}
	else{
		string te;
		switch(type){
			case LOAD_PER: te = "Loading: "; break;
			case INIT_PER: te = "Initialising: "; break;
			case RUN_PER: te = "Running: "; break;
			case RUN_GEN_PER: te = "Generation "+tstr(gen)+": "; break;
			case ANNEAL_PER: te = "Annealing: "; break;
			case OUTPUT_PER: te = "Outputting: "; break;
		}
		cout << te;
		cout.flush();
		
		/*
		string te;
		switch(type){
		case LOAD_PER: te = "Loading..."; break;
		case INIT_PER: te = "Initialising..."; break;
		case RUN_PER: te = "Running..."; break;
		case RUN_GEN_PER: te = "Generation "+tstr(gen)+"..."; break;
		case ANNEAL_PER: te = "Annealing..."; break;
		case OUTPUT_PER: te = "Outputting..."; break;
		}
		cout << te << endl;
		cout.flush();
		*/
		
		/*
		if(progress_on){ 
			progress << te << endl;
			progress.flush();
		}
		*/
	}
	
	percent_done = UNSET;
	percentage(0,1);
}


/// Ends the percentage
void percentage_end()
{
	percentage(1,1);
}


/// Prints percentage done
void percentage(double val, double val2)
{		
	if(!op()) return;
	
	if(!com_op && percent_done != LARGE){		
		if(percent_done == UNSET) percent_done = 0;
		
		/*
		// Draw progress bar like "Progress 0%...10%...20%..."
		auto per = 100*val/val2;
		while(percent_done <= per){
			if(int(percent_done)%10 == 0){
				cout << "" << percent_done << "%";
				if(percent_done == 100){ 
					percent_done = LARGE;
					cout << endl; cout.flush(); 
					return;
				}
			}
			cout << "."; 
			cout.flush();
			percent_done += PERCENT_STEP;
		}	
		*/

		auto per = 100*val/val2;
		while(percent_done <= per){
			if(int(percent_done)%10 == 0){
				printf("%d%%",(unsigned int)(percent_done));
				fflush(stdout); 
				//cout << percent_done << "%";
				
				/*
				if(progress_on){ 
					progress << percent_done << "%" << endl; 
					progress.flush();
				}
				*/
					
				if(percent_done == 100){ 
					percent_done = LARGE;
					cout << endl; cout.flush(); 
					return;
				}
			}
			printf(".");
			fflush(stdout); 
			//cout << "."; 
			cout.flush();
			percent_done += PERCENT_STEP;
		}		
		
		/*
		// Write progress values line-by-line to work around mpi's broken flush
		auto per = 100*val/val2;
		while(percent_done <= per){
			if(int(percent_done)%10 == 0){
				printf("%d%%",(unsigned int)(percent_done));
				fflush(stdout); 
				//cout << percent_done << "%";
						
				if(percent_done == 100){ 
					percent_done = LARGE;
					cout << endl; cout.flush(); 
					return;
				}
			}
			printf(".");
			fflush(stdout); 
			//cout << "."; 
			cout.flush();
			percent_done += PERCENT_STEP;
		}	
		*/		
	}
	else{
		auto per = (unsigned int)(100*val/val2);
		if(percent_done == UNSET || percent_done < per){
			cout << "<PROGRESS>" << per << endl;
			cout.flush();
			percent_done = per;
		}
	}
}


/// Finds the root for a file name
string find_root(string file)
{
	int i = file.length()-1; 
	while(i >= 0 && file.substr(i,1) != "\\" && file.substr(i,1) != "/") i--;

	return file.substr(0,i+1);
}


/// Waits indefinately
void wait()
{
	cout << "waiting" << endl; 
	while(true){};
}



/// Adds a fullstop at the end of a line
void add_full_stop(string &te)
{
	if(te.length() > 0){
		auto ch = te.substr(te.length()-1,1);
		if(ch != "." && ch != "?" && ch != "!" && ch != endli) te += ".";
	}
}


/// This read in a line to get information about memory usage
int parseLine(char* line)
{
	int i = strlen(line);
	const char* p = line;
	while (*p <'0' || *p > '9') p++;
	line[i-3] = '\0';
	i = atoi(p);
	return i;
}


/// Gets the memory being used by the process 
double memory_usage()
{
#ifdef WINDOWS
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;
	return virtualMemUsedByMe/1024.0;
#endif

#ifdef LINUX
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL){	
		if(strncmp(line, "VmSize:", 7) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
#endif

#ifdef MAC
	struct rusage usage;
	if(0 == getrusage(RUSAGE_SELF, &usage)){
		return usage.ru_maxrss/1024.0; // bytes
	}
	 else{
    		return 0;
	}
#endif

	return 0;
}


/// Returns total computer memory
double total_memory()
{
#ifdef WINDOWS
	MEMORYSTATUSEX memInfo;
	memInfo.dwLength = sizeof(MEMORYSTATUSEX);
	GlobalMemoryStatusEx(&memInfo);
	DWORDLONG totalVirtualMem = memInfo.ullTotalPageFile;
	return totalVirtualMem/1024.0;
#endif

#ifdef LINUX
	struct sysinfo memInfo;	
	sysinfo (&memInfo);
	long long totalVirtualMem = memInfo.totalram;
	totalVirtualMem += memInfo.totalswap;
	totalVirtualMem *= memInfo.mem_unit;
	totalVirtualMem /= 1024.0; // Convert to K
	return totalVirtualMem;
#endif

#ifdef MAC
	return LARGE;
#endif

	return LARGE;
}


/// The memory available 
double memory_available()
{
#ifdef WINDOWS
	MEMORYSTATUSEX memInfo;
	memInfo.dwLength = sizeof(MEMORYSTATUSEX);
	GlobalMemoryStatusEx(&memInfo);
	DWORDLONG availVirtualMem = memInfo.ullAvailPageFile;
	return availVirtualMem/1024.0;
#endif

#ifdef LINUX
	struct sysinfo memInfo;	
	sysinfo (&memInfo);
	return memInfo.freeram/1024.0;
#endif

#ifdef MAC
	return LARGE;
#endif

	return LARGE;
}


/// Outputs memory as a string
string mem_print(double mem)
{
	stringstream ss;
	ss << std::fixed;
  ss << setprecision(1);

	if(mem < 1024){   // Result
		ss << (unsigned int) mem << "K";
	}
	else{
		mem /= 1024;
		if(mem < 1024){ 
			if(mem < 10) ss << mem << "MB";
			else ss << (unsigned int) mem << "MB";
		}
		else{
			mem /= 1024;
			if(mem < 10) ss << mem << "GB";
			else ss << (unsigned int) mem << "GB";
		}
	}
	
	return ss.str();
}


/// Clips the exponential for individual effects (to ensure does go too large or small)
double exp_clip(double val)
{
	if(val < CLIP_MIN) val = CLIP_MIN;
	if(val > CLIP_MAX) val = CLIP_MAX;
	
	return exp(val);
}


/// Ensures the a value doesn't go below the threshold
double clip(double val, double thresh)
{
	if(val < thresh) return thresh;
	return val;
}


/// Gets ti (the division time) from an actual time
unsigned int get_ti(double t)
{
	return (unsigned int)(t);
}


/// Gets ti div time from div time (rounding down if equal)
unsigned int get_ti_lower(double t)
{
	auto ro = (unsigned int)(t);
	if(ro == t && ro != 0) return ro-1;
	return ro;
}


/// Gets next ti 
unsigned int get_ti_next(double t_next, const Details &details)
{
	//auto ti_next = (unsigned int)(ALMOST_ONE+((t_next-t_start)/dt));
	//if(ti_next > T) ti_next = T;
	auto ti_next = (unsigned int)(t_next);
	if(ti_next != t_next && ti_next < details.T) ti_next++;		
	
	return ti_next;
}


/// Gets ti, but ensures boundary
unsigned int get_ti_over(double tdiv)
{
	return (unsigned int)(tdiv+TINY);
}


/// Determines if time is on a timestep
bool on_timestep(double t, const Details &details) 
{
	double tdiv = calc_tdiv(t,details);
	auto ti = round_int(tdiv+TINY);
	
	if(dif(tdiv,ti,DIF_THRESH)) return false;
	return true;
}


/// Calculates tdiv from t
double calc_tdiv(double t, const Details &details)
{
	return (t - details.t_start)/details.dt;
}


/// Calculates t from tdiv
double calc_t(double tdiv, const Details &details)
{
	return details.t_start + tdiv*details.dt;
}


/// Adds algorithm warning to a list
void add_alg_warning(string te, unsigned int sample, vector <AlgWarn> &alg_warn)
{
	auto i = 0u; while(i < alg_warn.size() && alg_warn[i].te != te) i++;
	if(i < alg_warn.size()){
		alg_warn[i].num++;
	}
	else{
		AlgWarn aw;
		aw.te = te; aw.core = core(); aw.sample = sample; aw.num = 1;
		alg_warn.push_back(aw);
	}
}



/// Gets a string output from an equation type
void eqn_type_error(EqItemType type, unsigned int ref)
{
	string st;
	switch(type){
	case PARAMETER: st = "PARAMETER"; break;
	case PARAMVEC: st = "PARAMVEC"; break;
	case SPLINE: st = "SPLINE"; break;
	case SPLINEREF: st = "SPLINEREF"; break;
	case DERIVE: st = "DERIVE"; break;
	case IE: st = "IE"; break;
	case FE: st = "FE"; break;
	case POPNUM: st = "POPNUM"; break;
	case TIME: st = "TIME"; break;
	default: st = "DEFAULT"; break;
	}
	
	emsg("EQN TYPE ERROR: "+st+" "+tstr(ref));
}


/// Performs an operation on a set of numbers 
double calculate_operation(EqItemType op, vector <double> &num)
{
	auto N = num.size();
	
	switch(op){
	case ADD:
		{
			auto ans = 0.0; 
			for(auto val : num) ans += val;
			return ans;
		}
		
	case TAKE:
		if(N != 2){
			emsg("For TAKE should be 2");
		}
		return num[0]-num[1]; 
	
	case MULTIPLY:
		{
			auto ans = 1.0; 
			for(auto val : num) ans *= val;
			return ans;
		}
		
	case DIVIDE:
		if(N != 2) emsg("For DIVIDE should be 2");
		if(num[1] == 0) run_error("Equation caused a division by zero."); 
		return num[0]/num[1]; 
		
	case EXPFUNC:
		if(N != 1) emsg("For EXPFUNC should be 1");
		return exp(num[0]); 
			
	case SINFUNC:
		if(N != 1) emsg("For SINFUNC should be 1");
		return sin(num[0]); 
	
	case COSFUNC: 
		if(N != 1) emsg("For COSFUNC should be 1");
		return cos(num[0]); 
		
	case LOGFUNC: 
		if(N != 1) emsg("For LOGFUNC should be 1");
		if(num[0] <= 0) run_error("The quantity inside a log function became negative."); 
		return log(num[0]);
		
	case STEPFUNC:
		if(N != 1) emsg("For STEPFUNC should be 1");
		if(num[0] > 0) return 1;
		return 0; 
		
	case POWERFUNC:
		if(N != 2) emsg("For POWERFUNC should be 2");
		return pow(num[0],num[1]);
	
	case THRESHFUNC:
		if(N != 2) emsg("For THRESHFUNC should be 2");
		if(num[0] < num[1]) return 0;
		return num[0];
		
	case UBOUNDFUNC:
		if(N != 2) emsg("For UBOUNDFUNC should be 2");
		if(num[0] > num[1]) return INFINITY;
		return num[0];
	
	case MAXFUNC:
		if(N != 2) emsg("For MAXFUNC should be 2");
		if(num[0] > num[1]) return num[0];
		return num[1]; 
		
	case MINFUNC:
		if(N != 2) emsg("For MINFUNC should be 2");
		if(num[0] < num[1]) return num[0];
		return num[1]; 
		
	case ABSFUNC:
		if(N != 1) emsg("For MINFUNC should be 1");
		if(num[0] > 0) return num[0]; 
		return -num[0];
		
	case SQRTFUNC:
		if(N != 1) emsg("For SQRTFUNC should be 1");
		if(num[0] < 0) run_error("A square root of a negative number was found in an equation."); 
		return sqrt(num[0]);
		
	case SIGFUNC:
		if(N != 1) emsg("For SIGFUNCT should be 1");
		return 1/(1+exp(-num[0]));
		
	case NOOP: 
		return 0;
	
	default: emsg("Equation error7"); break;
	}
	
	return 0;
}


/// Tests the cubic spline algorithm 
void solve_cubic_spline_test()
{
	vector <double> x, f;
	 
	x.push_back(0); f.push_back(1);
	x.push_back(1); f.push_back(1.5);
	x.push_back(2); f.push_back(1.3);
	x.push_back(3); f.push_back(1.4);
	
	auto cspl = solve_cubic_spline(x,f,CUBIC_SPL);
	
	ofstream fout("../OP/spl.txt");
	
	for(auto xf = 0.0; xf <= 3; xf += 0.01){
		fout << xf << " " << calculate_cubic_spline(xf,cspl) << endl; 
	}
}


/// This provides an algorithm for solving cubic splines
CubicSpline solve_cubic_spline(const vector <double> &x, vector <double> f, SplineType type)
{
	CubicSpline cspl;
	cspl.type = type;
	
	auto n = x.size()-1;
	
	if(type == CUBICPOS_SPL){
		for(auto i = 0u; i <= n; i++){
			if(f[i] <= 0){
				run_error("For a 'Cubic +ve' spline all knot values must be positive");
			}
			
			f[i] = log(f[i]);
		}
	}
	
	vector <double> h(n), al(n), I(n), mu(n), z(n), a(n), b(n), c(n), d(n);
	
	for(auto i = 0u; i < n; i++){
		a[i] = f[i];
		h[i] = x[i+1]-x[i];
	}
	
	for(auto i = 1u; i < n; i++){
		al[i] = (3/h[i])*(f[i+1] - f[i]) - (3/h[i-1])*(f[i] - f[i-1]);
	}
	
	I[0] = 1; mu[0] = 0; z[0] = 0;
	
	for(auto i = 1u; i < n; i++){
		I[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1];
		mu[i] = h[i]/I[i];
		z[i] = (al[i]-h[i-1]*z[i-1])/I[i];
	}
	
	for(int j = n-1; j >= 0; j--){
		auto cne = 0.0; if(j < int(n-1)) cne = c[j+1];
		c[j] = z[j] - mu[j]*cne;
		b[j] = ((f[j+1]-f[j])/h[j]) - ((h[j]*(cne+2*c[j])/3));
		d[j] = (cne-c[j])/(3*h[j]);
	}
	
	cspl.x = x; cspl.a = a; cspl.b = b; cspl.c = c; cspl.d = d;
	
	return cspl;
}


/// Calculates a cubic spline from 
double calculate_cubic_spline(double xf, const CubicSpline &cspl) 
{
	const auto &x = cspl.x;
	 
	if(xf < x[0]) run_error("Value '"+tstr(xf)+"' outside of spline range");
	
	auto i = 0u;
	while(i+1 < x.size() && x[i+1] < xf) i++;
	
	if(i+1 == x.size()) run_error("Value '"+tstr(xf)+"' outside of spline range");
	
	auto dx = xf-x[i];
	
	auto val = cspl.a[i] + cspl.b[i]*dx + cspl.c[i]*dx*dx + cspl.d[i]*dx*dx*dx;
	
	if(cspl.type == CUBICPOS_SPL) return exp(val);
	return val;
}


/// Calculates a cubic spline from 
void calculate_cubic_spline_precalc(vector <double> &precalc, unsigned int i, const vector <unsigned int> &list_time, const vector <CubicDiv> &div, const CubicSpline &cspl)
{
	const auto &a = cspl.a, &b = cspl.b, &c = cspl.c, &d = cspl.d;
	
	if(cspl.type == CUBICPOS_SPL){ 
		for(auto ti : list_time){
			const auto &di = div[ti];
			auto k = di.i;
			precalc[i+ti] = exp(a[k] + b[k]*di.dx + c[k]*di.dx2 + d[k]*di.dx3);
		}
	}
	else{
		for(auto ti : list_time){
			const auto &di = div[ti];
			auto k = di.i;
			precalc[i+ti] = a[k] + b[k]*di.dx + c[k]*di.dx2 + d[k]*di.dx3;
		
		}
	}
}
