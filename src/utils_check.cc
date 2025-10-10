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

using namespace std;

#include "const.hh"
#include "utils.hh"
#include "matrix.hh"

/// Calculates the tri_gamma function
double tri_gamma(double z)
{
	auto sum = 0.0;
	//for(auto n = 0u; n < LARGE; n++){
	//for(auto n = 0u; n < 1000*z; n++){
	for(auto n = 0u; n < 10000000; n++){
			
		auto d = 1.0/((z+n)*(z+n));
		sum += d;
	}
	
	/*
	auto dx = 0.01;
	auto sum2 = 0.0;
	for(auto x = dx; x < 1; x+= dx){
		sum2 += dx*pow(x,z-1)*log(x)/(1-x);
	}
	*/
	
	//return -sum2;
	//cout << sum << " " << sum2 << " j" << endl;
	return sum;
}


/// Calculates the tri_gamma function
double di_gamma(double z)
{
	auto c = 0.577215664901532860;
	
	auto sum = -c;
	//for(auto n = 0u; n < LARGE; n++){
	for(auto n = 0u; n < 10*z; n++){
		auto d = (z-1)/((1+n)*(z+n));
		sum += d;
	}	
	
	return sum;
}


/// This tests if the Jeffreys expression in the paper are correct or not
void test_jeffreys()
{
	if(false){  // exp_mean
		auto pmin = 0.05, pmax = 5.0, dp = 0.01;
		auto xmax = 10000.0, dx = 0.01;
		
		ofstream fout("Output/mean.txt");
		for(auto p = pmin; p < pmax; p += 0.1){
			auto sum = 0.0, sum2 = 0.0;
			for(auto x = dx; x < xmax; x += dx){
				auto v_pp = exp_mean_probability(x,p+dp);
				auto v_pm = exp_mean_probability(x,p-dp);
				auto v = exp_mean_probability(x,p);
				auto dpp = (v_pp+v_pm-2*v)/(dp*dp);
				auto f = exp(v);
				sum += dpp*f*dx;
				sum2 += f*dx;
			}
			sum /= sum2;
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			if(sum < 0) sum = -sum;
			fout << p << " " << sqrt(sum) << " " << 0.5/p << endl; 
		}
	}
	
	if(false){  // exp_rate
		auto pmin = 0.05, pmax = 5.0, dp = 0.01;
		auto xmax = 10000.0, dx = 0.01;
		
		ofstream fout("Output/rate.txt");
		for(auto p = pmin; p < pmax; p += 0.1){
			auto sum = 0.0, sum2 = 0.0;
			for(auto x = dx; x < xmax; x += dx){
				auto v_pp = exp_rate_probability(x,p+dp);
				auto v_pm = exp_rate_probability(x,p-dp);
				auto v = exp_rate_probability(x,p);
				auto dpp = (v_pp+v_pm-2*v)/(dp*dp);
				auto f = exp(v);
				sum += dpp*f*dx;
				sum2 += f*dx;
			}
			sum /= sum2;
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			if(sum < 0) sum = -sum;
			fout << p << " " << sqrt(sum) << " " << 0.5/p << endl; 
		}
	}
	
	if(false){  // gamma fixed
		auto m = 1.0;
		auto pmin = 0.05, pmax = 2.0, dp = 0.001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/gamma_fixed.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum = 0.0, sum2 = 0.0;
			for(auto x = 0.1*dx; x < xmax; x += dx){
				auto v_pp = gamma_probability(x,m,p+dp);
				auto v_pm = gamma_probability(x,m,p-dp);
				auto v = gamma_probability(x,m,p);
				auto dpp = (v_pp+v_pm-2*v)/(dp*dp);
				auto f = exp(v);
				sum += dpp*f*dx;
				sum2 += f*dx;
			}
			sum /= sum2;
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			if(sum < 0) sum = -sum;
			fout << p << " " << sqrt(sum) << " " << sqrt(tri_gamma(1/(p*p)) - (p*p))/(p*p*p) << " " << 1/(p) << endl; 
		}
	}
	
	if(false){  // gamma 
		auto m = 1.0, dm = 0.001;
		auto pmin = 0.05, pmax = 2.0, dp = 0.001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/gamma.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum2 = 0.0;
			auto dmm_sum = 0.0, dmp_sum = 0.0, dpp_sum = 0.0;
			for(auto x = 0.1*dx; x < xmax; x += dx){
				auto v_pp = gamma_probability(x,m+dm,p+dp);
				auto v_pm = gamma_probability(x,m+dm,p-dp);
				auto v_mp = gamma_probability(x,m-dm,p+dp);
				auto v_mm = gamma_probability(x,m-dm,p-dp);
				auto v = gamma_probability(x,m,p);
				auto dmm = ((v_pp+v_pm)/2 + (v_mp+v_mm)/2 - 2*v)/(dm*dm);
				auto dmp = (v_pp+v_mm - v_pm - v_mp)/(4*dm*dp);
				auto dpp = ((v_pp+v_mp)/2 + (v_mm+v_pm)/2 - 2*v)/(dp*dp);
				auto f = exp(v);
				dmm_sum += dmm*f*dx;
				dmp_sum += dmp*f*dx;
				dpp_sum += dpp*f*dx;
				sum2 += f*dx;
			}
			dmm_sum /= sum2;
			dmp_sum /= sum2;
			dpp_sum /= sum2;
			
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			auto det = dmm_sum*dpp_sum - dmp_sum*dmp_sum;
				
			//cout << tri_gamma(1/(p*p)) << " " << (p*p) << "U" << endl;
			fout << p << " " << sqrt(det) << " " << sqrt(tri_gamma(1/(p*p)) - (p*p))/(p*p*p*p) << " " << 1/(p*p) << endl; 
		}
	}
	
	
	if(false){  // lognormal fixed
		auto m = 1.0;
		auto pmin = 0.05, pmax = 2.0, dp = 0.001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/lognormal_fixed.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum = 0.0, sum2 = 0.0;
			for(auto x = 0.1*dx; x < xmax; x += dx){
				auto v_pp = lognormal_probability(x,m,p+dp);
				auto v_pm = lognormal_probability(x,m,p-dp);
				auto v = lognormal_probability(x,m,p);
				auto dpp = (v_pp+v_pm-2*v)/(dp*dp);
				auto f = exp(v);
				sum += dpp*f*dx;
				sum2 += f*dx;
			}
			sum /= sum2;
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			if(sum < 0) sum = -sum;
			fout << p << " " << sqrt(sum) << " " << p/(log(1+p*p)*(1+p*p)) << " " << 2/(p) << endl; 
		}
	}
	
	if(false){  // lognormal 
		auto m = 1.0, dm = 0.0001;
		auto pmin = 0.05, pmax = 2.0, dp = 0.0001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/lognormal.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum2 = 0.0;
			auto dmm_sum = 0.0, dmp_sum = 0.0, dpp_sum = 0.0;
			for(auto x = 0.1*dx; x < xmax; x += dx){
				auto v_pp = lognormal_probability(x,m+dm,p+dp);
				auto v_pm = lognormal_probability(x,m+dm,p-dp);
				auto v_mp = lognormal_probability(x,m-dm,p+dp);
				auto v_mm = lognormal_probability(x,m-dm,p-dp);
				auto v = lognormal_probability(x,m,p);
				auto dmm = ((v_pp+v_pm)/2 + (v_mp+v_mm)/2 - 2*v)/(dm*dm);
				auto dmp = (v_pp+v_mm - v_pm - v_mp)/(4*dm*dp);
				auto dpp = ((v_pp+v_mp)/2 + (v_mm+v_pm)/2 - 2*v)/(dp*dp);
				auto f = exp(v);
				dmm_sum += dmm*f*dx;
				dmp_sum += dmp*f*dx;
				dpp_sum += dpp*f*dx;
				sum2 += f*dx;
			}
			dmm_sum /= sum2;
			dmp_sum /= sum2;
			dpp_sum /= sum2;
			
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			auto det = dmm_sum*dpp_sum - dmp_sum*dmp_sum;
				
			//cout << tri_gamma(1/(p*p)) << " " << (p*p) << "U" << endl;
			fout << p << " " << sqrt(det) << " " << p/(pow(log(1+p*p),1.5)*(1+p*p)) << " " << 2/(p*p) << endl; 
		}
	}
	
	
	if(false){  // weibull fixed
		auto m = 1.0;
		//auto pmin = 0.05, pmax = 2.0, dp = 0.001;
		auto pmin = 0.5, pmax = 40.0, dp = 0.001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/weibull_fixed.txt");
		for(auto p = pmin; p < pmax; p += 0.1){
			auto sum = 0.0, sum2 = 0.0;
			for(auto x = 0.1*dx; x < xmax; x += dx){
				auto v_pp = weibull_probability(x,m,p+dp);
				auto v_pm = weibull_probability(x,m,p-dp);
				auto v = weibull_probability(x,m,p);
				auto dpp = (v_pp+v_pm-2*v)/(dp*dp);
				auto f = exp(v);
				sum += dpp*f*dx;
				sum2 += f*dx;
			}
			sum /= sum2;
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			if(sum < 0) sum = -sum;
			fout << p << " " << sqrt(sum) << " " << p/(log(1+p*p)*(1+p*p)) << " " << 2/(p) << endl; 
		}
	}
	
	if(false){  // weibull
		auto m = 1.0, dm = 0.001;
		//auto pmin = 0.05, pmax = 2.0, dp = 0.001;
		auto pmin = 0.5, pmax = 40.0, dp = 0.001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/weibull.txt");
		for(auto p = pmin; p < pmax; p += 0.1){
			auto sum2 = 0.0;
			auto dmm_sum = 0.0, dmp_sum = 0.0, dpp_sum = 0.0;
			for(auto x = 0.1*dx; x < xmax; x += dx){
				auto v_pp = weibull_probability(x,m+dm,p+dp);
				auto v_pm = weibull_probability(x,m+dm,p-dp);
				auto v_mp = weibull_probability(x,m-dm,p+dp);
				auto v_mm = weibull_probability(x,m-dm,p-dp);
				auto v = weibull_probability(x,m,p);
				auto dmm = ((v_pp+v_pm)/2 + (v_mp+v_mm)/2 - 2*v)/(dm*dm);
				//cout << v_pp << " " << v_pm << " " << v_mp  << " "<< v_mm << "uu" << endl;
				auto dmp = (v_pp+v_mm - v_pm - v_mp)/(4*dm*dp);
				auto dpp = ((v_pp+v_mp)/2 + (v_mm+v_pm)/2 - 2*v)/(dp*dp);
				auto f = exp(v);
				//cout << f << "f" << endl;
				dmm_sum += dmm*f*dx;
				dmp_sum += dmp*f*dx;
				dpp_sum += dpp*f*dx;
				sum2 += f*dx;
			}
			dmm_sum /= sum2;
			dmp_sum /= sum2;
			dpp_sum /= sum2;
			
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			auto det = dmm_sum*dpp_sum - dmp_sum*dmp_sum;
				
			//cout << tri_gamma(1/(p*p)) << " " << (p*p) << "U" << endl;
			fout << p << " " << sqrt(det) << " " << p/(pow(log(1+p*p),1.5)*(1+p*p)) << " " << 2/(p*p) << endl; 
		}
	}
	
	
	if(false){  // normal fixed
		auto m = 1.0;
		auto pmin = 0.05, pmax = 2.0, dp = 0.001;
		auto xmax = 100.0, dx = 0.0001;
		
		ofstream fout("Output/normal_fixed.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum = 0.0, sum2 = 0.0;
			for(auto x = -xmax; x < xmax; x += dx){
				auto v_pp = normal_probability(x,m,p+dp);
				auto v_pm = normal_probability(x,m,p-dp);
				auto v = normal_probability(x,m,p);
				auto dpp = (v_pp+v_pm-2*v)/(dp*dp);
				auto f = exp(v);
				sum += dpp*f*dx;
				sum2 += f*dx;
			}
			sum /= sum2;
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			if(sum < 0) sum = -sum;
			fout << p << " " << sqrt(sum) << " " << 2/p << " " << 2/(p) << endl; 
		}
	}
	
	if(false){  // normal 
		auto m = 1.0, dm = 0.0001;
		auto pmin = 0.05, pmax = 2.0, dp = 0.0001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/normal.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum2 = 0.0;
			auto dmm_sum = 0.0, dmp_sum = 0.0, dpp_sum = 0.0;
			for(auto x = -xmax; x < xmax; x += dx){
				auto v_pp = normal_probability(x,m+dm,p+dp);
				auto v_pm = normal_probability(x,m+dm,p-dp);
				auto v_mp = normal_probability(x,m-dm,p+dp);
				auto v_mm = normal_probability(x,m-dm,p-dp);
				auto v = normal_probability(x,m,p);
				auto dmm = ((v_pp+v_pm)/2 + (v_mp+v_mm)/2 - 2*v)/(dm*dm);
				auto dmp = (v_pp+v_mm - v_pm - v_mp)/(4*dm*dp);
				auto dpp = ((v_pp+v_mp)/2 + (v_mm+v_pm)/2 - 2*v)/(dp*dp);
				auto f = exp(v);
				dmm_sum += dmm*f*dx;
				dmp_sum += dmp*f*dx;
				dpp_sum += dpp*f*dx;
				sum2 += f*dx;
			}
			dmm_sum /= sum2;
			dmp_sum /= sum2;
			dpp_sum /= sum2;
			
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			auto det = dmm_sum*dpp_sum - dmp_sum*dmp_sum;
				
			//cout << tri_gamma(1/(p*p)) << " " << (p*p) << "U" << endl;
			fout << p << " " << sqrt(det) << " " << 2/(p*p) << " " << 2/(p*p) << endl; 
		}
	}
	
	/*
	if(false){  // mvn fixed
		auto m = 1.0;
		auto pmin = 0.05, pmax = 2.0, dp = 0.001;
		auto xmax = 100.0, dx = 0.0001;
		vector < vector <double> > covar;
		covar.resize(2);
		covar.push_back(1);	covar.push_back(0);
		covar.push_back(0);	covar.push_back(1);
		
		vector <double> mean;
		mean.push_back(0); mean.push_back(0); 
		
		
		ofstream fout("Output/mvn_fixed.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum = 0.0, sum2 = 0.0;
			for(auto xi = -xmax; xi < xmax; xi += dx){
				for(auto xj = -xmax; xj < xmax; xj += dx){
					vector <double> x;
					x.push_back(xi); x.push_back(xj);
					
				auto v_pp = mvn_probability(x,m,p+dp);
				auto v_pm = normal_probability(x,m,p-dp);
				auto v = normal_probability(x,m,p);
				auto dpp = (v_pp+v_pm-2*v)/(dp*dp);
				auto f = exp(v);
				sum += dpp*f*dx;
				sum2 += f*dx;
			}
			sum /= sum2;
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			if(sum < 0) sum = -sum;
			fout << p << " " << sqrt(sum) << " " << 2/p << " " << 2/(p) << endl; 
		}
	}
	
	if(false){  // normal 
		auto m = 1.0, dm = 0.0001;
		auto pmin = 0.05, pmax = 2.0, dp = 0.0001;
		auto xmax = 100.0, dx = 0.001;
		
		ofstream fout("Output/normal.txt");
		for(auto p = pmin; p < pmax; p += 0.01){
			auto sum2 = 0.0;
			auto dmm_sum = 0.0, dmp_sum = 0.0, dpp_sum = 0.0;
			for(auto x = -xmax; x < xmax; x += dx){
				auto v_pp = normal_probability(x,m+dm,p+dp);
				auto v_pm = normal_probability(x,m+dm,p-dp);
				auto v_mp = normal_probability(x,m-dm,p+dp);
				auto v_mm = normal_probability(x,m-dm,p-dp);
				auto v = normal_probability(x,m,p);
				auto dmm = ((v_pp+v_pm)/2 + (v_mp+v_mm)/2 - 2*v)/(dm*dm);
				auto dmp = (v_pp+v_mm - v_pm - v_mp)/(4*dm*dp);
				auto dpp = ((v_pp+v_mp)/2 + (v_mm+v_pm)/2 - 2*v)/(dp*dp);
				auto f = exp(v);
				dmm_sum += dmm*f*dx;
				dmp_sum += dmp*f*dx;
				dpp_sum += dpp*f*dx;
				sum2 += f*dx;
			}
			dmm_sum /= sum2;
			dmp_sum /= sum2;
			dpp_sum /= sum2;
			
			if(sum2 < 0.99) cout << p << " " << sum2 <<  " ran" << endl;
			
			auto det = dmm_sum*dpp_sum - dmp_sum*dmp_sum;
				
			//cout << tri_gamma(1/(p*p)) << " " << (p*p) << "U" << endl;
			fout << p << " " << sqrt(det) << " " << 2/(p*p) << " " << 2/(p*p) << endl; 
		}
	}
	*/
	
	/*
	{
		ofstream fout("Output/trigamma.txt");
		for(auto cv = 0.05; cv < 1; cv += 0.1){
			fout << cv << " " << tri_gamma(1/(cv*cv)) << " " << cv*cv <<  endl;
		}
	}
	
	{
		ofstream fout("Output/digamma.txt");
		for(auto cv = 0.05; cv < 1; cv += 0.1){
			fout << cv << " " << di_gamma(1/(cv*cv)) << " " << cv*cv <<  endl;
		}
	}
	*/
	cout << "done" << endl;
}

	
/// Checks that a distribution is being generated correctly
void test_distribution()
{
	string warn;
	
	{
		for(auto shape = 1.0; shape < 5; shape += 0.1){
			auto gam1 = exp(lgamma(1+2.0/shape));
			auto gam2 = exp(lgamma(1+1.0/shape));
			
			auto cv = sqrt(gam1 - gam2*gam2)/gam2;
			cout << shape << " "<< cv << " " << 1/sqrt(shape) << " hh" << endl;
		}
	
		return;
	}
	
	{
		auto lam = 10;
		
		auto sum = 0.0;
		for(auto i = 0u; i < 100; i++){
			sum += exp(poisson_probability(i,lam));
			cout << 1-sum << " " << poisson_upper_probability_no_log(i,lam) << "comp" << endl;
		}
	}
	return;
		
	{
		auto m = 10.0, p = 0.5;
		auto sum = 0.0;
		for(auto k = 0u; k < 100; k++){
			auto val = neg_binomial_probability(k,m,p);
			cout << k << " " << val << endl;
			sum += exp(val);
		}			
		cout << sum << " sum" << endl;
		emsg("check dist");
	}
	
	{
		auto m = 30.0, cv = 0.5;
		
		auto dx = 0.0001;
		auto sum = 0.0;
		
		auto avd = 0.0, avd2 = 0.0;
		for(auto x = 0.0; x < 100; x += dx){
			sum += exp(gamma_probability(x,m,cv))*dx;
			if(ran() < 0.001){
				cout << 1-sum << " " << gamma_upper_probability_no_log(x,m,cv) << " " << exp(gamma_upper_probability(x,m,cv) ) << " comp" << endl;
			}
			
			avd += x*exp(gamma_probability(x,m,cv))*dx;
			avd2 += x*x*exp(gamma_probability(x,m,cv))*dx;
		}
		cout << sum << "sum" << endl;
		
		auto av = 0.0, av2 = 0.0;
		auto loopmax = 1000000u;
		for(auto loop = 0u; loop < loopmax; loop++){
			auto x = gamma_sample(m,cv,warn);
			av += x;
			av2 += x*x;
		}
		av /= loopmax; av2 /= loopmax;
		auto shape = 1.0/(cv*cv);
		cout << 30*30/shape << endl;
		cout << avd << " " << avd2-avd*avd <<" h12" << endl;
		cout << av << " " << av2-av*av <<" h12" << endl;
	}
	
	return;
	
	auto imax = 100000u;
	auto sum = 0.0, sum2 = 0.0; 
	auto meansamp = 3.0, cvsamp = 1.5, sdsamp=2.0;
	auto shapesamp = 2.0, scalesamp=3.0, lamsamp=5.0, zsamp = 0.5;
	auto alphasamp = 0.5, betasamp = 2.0;
	
	auto dx = 0.00001;	
	auto dsum = 0.0, dsum2 = 0.0, ds = 0.0;

	switch(6){
	case 0: // Tests the normal distribution
		{
			for(auto i = 0u; i < imax; i++){
				auto val = normal_sample(meansamp,sdsamp,warn);
				sum += val; sum2 += val*val;
			}

			for(auto x = -10.0; x < 10; x+=dx){
				auto val = exp(normal_probability(x,meansamp,sdsamp));
				dsum += val*x; dsum2 += val*x*x; ds += val;
			}
		
			auto mean = sum/imax, var = sum2/imax - mean*mean;
			
			cout << meansamp << " "<< sdsamp*sdsamp << "True" << endl;
			cout << mean << " " << var << " Sample"  << endl;
			cout << dsum/ds << " " << (dsum2/ds - (dsum/ds)*(dsum/ds)) << " Distribution"  << endl;
		}
		break;
	
	case 1: // Tests the gamma distribution
		{
			for(auto i = 0u; i < imax; i++){
				auto val = gamma_sample(meansamp,cvsamp,warn);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 100; x+=dx){
				auto val = exp(gamma_probability(x,meansamp,cvsamp));
				dsum += val*x; dsum2 += val*x*x; ds += val;
			}
	
			auto mean = sum/imax, var = sum2/imax - mean*mean;
			auto sd = sqrt(var), cv = sd/mean;
	
			cout << meansamp << " "<< cvsamp << "True gamma" << endl;
			cout << mean << " " << cv << " Sample"  << endl;
			
			auto var_dist = dsum2/ds - (dsum/ds)*(dsum/ds);
			auto mean_dist = dsum/ds;
			auto cv_dist = sqrt(var_dist)/mean_dist;
			
			cout << mean_dist << " " << cv_dist << " Distribution"  << endl;
		}
		break;
		
	case 2: // Tests the log-normal distribution
		{
			for(auto i = 0u; i < imax; i++){
				auto val = lognormal_sample(meansamp,cvsamp,warn);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 100; x+=dx){
				auto val = exp(lognormal_probability(x,meansamp,cvsamp));
				dsum += val*x; dsum2 += val*x*x; ds += val;
			}
	
			auto mean = sum/imax, var = sum2/imax - mean*mean;
			auto sd = sqrt(var), cv = sd/mean;
	
			cout << meansamp << " "<< cvsamp << "True lognormal" << endl;
			cout << mean << " " << cv << " Sample"  << endl;
			
			auto var_dist = dsum2/ds - (dsum/ds)*(dsum/ds);
			auto mean_dist = dsum/ds;
			auto cv_dist = sqrt(var_dist)/mean_dist;
			
			cout << mean_dist << " " << cv_dist << " Distribution"  << endl;
		}
		break;
		
	case 3: // Tests the Weibull distribution
		{
			for(auto i = 0u; i < imax; i++){
				auto val = weibull_sample(scalesamp,shapesamp,warn);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 100; x+=dx){
				auto val = exp(weibull_probability(x,scalesamp,shapesamp));
				dsum += val*x; dsum2 += val*x*x; ds += val;
			}
	
			auto mean = sum/imax, var = sum2/imax - mean*mean;
			
			auto meansa = scalesamp*tgamma(1+1.0/shapesamp);
			auto varsa = scalesamp*scalesamp*(tgamma(1+2.0/shapesamp) - tgamma(1+1.0/shapesamp)*tgamma(1+1.0/shapesamp));
			cout << meansa << " "<< varsa << "True weibull" << endl;
			cout << mean << " " << var << " Sample"  << endl;
			
			auto var_dist = dsum2/ds - (dsum/ds)*(dsum/ds);
			auto mean_dist = dsum/ds;
		
			cout << mean_dist << " " << var_dist << " Distribution"  << endl;
		}
		break;
		
	case 4: // Tests the Poisson distribution
		{
			for(auto i = 0u; i < imax; i++){
				auto val = poisson_sample(lamsamp,warn);
				sum += val; sum2 += val*val;
			}

			for(auto k = 0u; k < 100; k++){
				auto val = exp(poisson_probability(k,lamsamp));
				dsum += val*k; dsum2 += val*k*k; ds += val;
			}
	
			auto mean = sum/imax, var = sum2/imax - mean*mean;
			
			cout << lamsamp << " "<< lamsamp << "True Poisson" << endl;
			cout << mean << " " << var << " Sample"  << endl;
			
			auto var_dist = dsum2/ds - (dsum/ds)*(dsum/ds);
			auto mean_dist = dsum/ds;
		
			cout << mean_dist << " " << var_dist << " Distribution"  << endl;
		}
		break;
		
	case 5: // Tests the bernoulli distribution
		{
			for(auto i = 0u; i < imax; i++){
				auto val = bernoulli_sample(zsamp,warn);
				sum += val; sum2 += val*val;
			}

			for(auto x = 0u; x < 2; x++){
				auto val = exp(bernoulli_probability(x,zsamp));
				dsum += val*x; dsum2 += val*x*x; ds += val;
			}
	
			auto mean = sum/imax, var = sum2/imax - mean*mean;
		
			cout << zsamp << " " << zsamp*(1-zsamp) << "True bernoulli" << endl;
			cout << mean << " " << var << " Sample"  << endl;
			
			auto mean_dist = dsum/ds;
			auto var_dist = dsum2/ds - (dsum/ds)*(dsum/ds);
			
			cout << mean_dist << " " << var_dist << " Distribution"  << endl;
		}
		break;
		
	case 6: // Tests the beta distribution
		{
			for(auto i = 0u; i < imax; i++){
				auto val = beta_sample(alphasamp,betasamp,warn);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 1-dx; x+=dx){
				auto val = exp(beta_probability(x,alphasamp,betasamp,warn));
				dsum += val*x; dsum2 += val*x*x; ds += val;
			}
	
			auto mean = sum/imax, var = sum2/imax - mean*mean;
			auto meansa = alphasamp/( alphasamp+betasamp);
			auto ab =(alphasamp+betasamp);
			auto varsa = alphasamp*betasamp/(ab*ab*(ab+1));
			
			cout << meansa << " "<< varsa << "True beta" << endl;
			cout << mean << " " << var << " Sample"  << endl;
			
			auto var_dist = dsum2/ds - (dsum/ds)*(dsum/ds);
			auto mean_dist = dsum/ds;
			
			cout << mean_dist << " " << var_dist << " Distribution"  << endl;
		}
		break;
		
	}
}


/// Generates synthetic data
void generate_data()
{	
	ofstream fout("Examples/loc.csv");	
	{
		auto d = 11;
		ofstream fout("Examples/loc.csv");
		fout << "ID,x,y,col" << endl;
		for(auto j = 0u; j < 5; j++){
			for(auto i = 0u; i < 5; i++){
				fout << "L-" << (j+1) << "-" << i+1 << "," << i*d << "," << j*d << ",#cccccc" << endl;
			}
		}
	}
	
	{
		auto D = 50;
		ofstream fout("Examples/loc2.csv");
		fout << "ID,x,y,col" << endl;
		for(auto j = 0u; j < 20; j++){
			auto x = ran()*D, y = ran()*D;
			fout << "F-" << (j+1) << "," << x << "," << y << ",#55dd55" << endl;
		}
		for(auto j = 0u; j < 20; j++){
			auto x = ran()*D, y = ran()*D;
			fout << "S-" << (j+1) << "," << x << "," << y << ",#5555dd" << endl;
		}
	}
	
	{
		ofstream fout("Examples/ind.csv");
		fout << "ID,t,DS" << endl;
		for(auto i = 0u; i < 50; i++){
			fout << "Ind-" << i << ",0,S" << endl;
			fout << endl;
		}
	}
	
	{
		ofstream fout("Examples/move.csv");
		fout << "ID,t,from,to" << endl;
		fout << "Ind-0,S,I" << endl; 
	}
	
	{
		ofstream fout("Examples/cow.csv");
		fout << "ID,t,DS" << endl;
		for(auto i = 0u; i < 50; i++){
			fout << "Cow-" << i << ",0,Sc" << endl;
			fout << endl;
		}
	}
		
	{
		ofstream fout("Examples/bad.csv");
		fout << "ID,t,DS" << endl;
		for(auto i = 0u; i < 50; i++){
			fout << "Bad-" << i << ",0,Sb" << endl;
			fout << endl;
		}
	}
	
	{
		ofstream fout("Examples/bad-move.csv");
		fout << "ID,t,DS" << endl;
		fout << "Bad-0,Sb,Ib" << endl; 
	}
	return;
	
	{
		ofstream fout("Examples/individual_data.csv");
	
		string warn;
		fout << "ID,X,DS,Sex" << endl;
		for(auto i = 0u; i < 500; i++){
			fout << "Ind-" << i << "," << 10+normal_sample(0,1,warn) << ","; 
			if(i < 2) fout << "I"; else fout << "S";
			fout << ",";
			if(ran() < 0.5) fout << "M"; else fout << "F";
			fout << endl;
		}
		cout << "Data generated" << endl;
	}
}


/// Simulates a set of individuals for a disease transmission experiment
void simulate_trans_exp()
{
	auto npar = 50u;
	auto npro = 20u;
		
	auto N = npar*npro;
	auto G = 10;                            // Number of individual per group
		
	auto n = npar + N;
		
	vector <string> name, name_in;
	for(auto i = 0u; i < n; i++){
		string na = "Ind."+to_string(i);
		if(i < npar) na = "Par."+to_string(i);
		else na = "Ind."+to_string(i);
		
		name.push_back(na);
		if(i >= npar) name_in.push_back(na);
	}
	
	ofstream fout2("Examples/ind.csv");
	fout2 << "ID,DS,Group" << endl;

	for(auto i = npar; i < n; i++){
		cout << i << " " << name_in.size() << "ai" << endl;
		auto j = (unsigned int)(name_in.size()*ran());
		
		fout2 << name_in[j];
		name_in.erase(name_in.begin()+j);
		
		auto z = (i-npar)/G;
	
		auto k = i-npar-z*G;
		if(k < 2) fout2 << ",E"; 
		else fout2 << ",S";
		
		fout2 << "," << "Gr" << z+1;
		fout2 << endl;
	}

	ofstream fout("Examples/A.csv");
	for(auto i = 0u; i < n; i++){
		if(i != 0) fout << ",";
		fout << name[i];
	}
	fout << endl;
	
	for(auto j = 0u; j < n; j++){
		for(auto i = 0u; i < n; i++){
			if(i != 0) fout << ",";
			if(i == j) fout << "1";
			else{
				if(j < npar){
					if(i < npar){
						fout << "0";
					}
					else{
						auto p = (i-npar)/npro;
						if(j == p) fout << "0.5"; else fout << "0";
					}
				}
				else{
					if(i < npar){
						auto p = (j-npar)/npro;
						if(i == p) fout << "0.5"; else fout << "0";
					}
					else{
						auto p1 = (j-npar)/npro;
						auto p2 = (i-npar)/npro;
						if(p1 == p2) fout << "0.25"; else fout << "0";
					}
				}
			}
		}
		fout << endl;
	}
	
	{
		string warn;
		
		ofstream fout("Examples/X.csv");
		fout << "Ind,Value" << endl;
		for(auto i = 0u; i < n; i++){
			fout << name[i] << "," << normal_sample(0,1,warn) << endl;
		}
	}
}


/// Tests to see if incomplete distributions work
void test_incomplete_distribution()
{ 
	auto mean = 2.0, cv = 0.3;
	
	{
		cout << "start" << endl;
		auto sum = 0.0;
		for(auto loop = 0u; loop < 100000000; loop++){
			//auto val = gamma_probability(1,mean,cv);
			//sum += val;
			
			//for(auto k = 0u; k <  150; k++) sum += 1;
		
			for(auto k = 0u; k < 4; k++){
				//sum += 1;
				sum += log(1+k);
			}
		}
		cout << sum << "end" << endl;
		return;
	}
	
	auto dx = 0.0001;
	auto sum = 0.0;
	vector <double> sum_st;
	
	switch(4){
	case 0: // Test gamma distribution
		{
			for(auto x = 0.0; x < 10; x += dx){
				auto val = exp(gamma_probability(x,mean,cv));
				sum_st.push_back(sum);
				sum += val*dx;
			}
			
			auto k = 0u;
			for(auto x = 0.0; x < 10; x += dx){
				if(k%1000 == 0){
					cout << x << " " << log(1-sum_st[k]/sum) << " " <<  gamma_upper_probability(x,mean,cv) << endl;
				}
				k++;
			}	
		}
		break;
		
	case 1: // Tests log-normal distribution
		{
			for(auto x = 0.0; x < 10; x += dx){
				auto val = exp(lognormal_probability(x,mean,cv));
				sum_st.push_back(sum);
				sum += val*dx;
			}
			
			auto k = 0u;
			for(auto x = 0.0; x < 10; x += dx){
				if(k%1000 == 0){
					cout << x << " " << log(1-sum_st[k]/sum) << " " << lognormal_upper_probability(x,mean,cv) << endl;
				}
				k++;
			}	
		}
		break;
		
	case 2: // Tests Weibull distribution
		{
			auto scale = 1.5;
			auto shape = 2;
			
			for(auto x = 0.0; x < 10; x += dx){
				auto val = exp(weibull_probability(x,scale,shape));
				sum_st.push_back(sum);
				sum += val*dx;
			}
			
			auto k = 0u;
			for(auto x = 0.0; x < 10; x += dx){
				if(k%1000 == 0){
					cout << x << " " << log(1-sum_st[k]/sum) << " " << weibull_upper_probability(x,scale,shape) << endl;
				}
				k++;
			}	
		}
		break;
		
	case 3: // Tests exp rate distribution
		{
			//auto scale = 1.5;
			auto rate = 0.5;
			
			for(auto x = 0.0; x < 100; x += dx){
				auto val = exp(exp_rate_probability(x,rate));
				sum_st.push_back(sum);
				sum += val*dx;
			}
			
			auto k = 0u;
			for(auto x = 0.0; x < 100; x += dx){
				if(k%1000 == 0){
					cout << x << " " << log(1-sum_st[k]/sum) << " " << exp_rate_upper_probability(x,rate) << endl;
				}
				k++;
			}	
		}
		break;
		
	case 4: // Tests exp mean distribution
		{	
			for(auto x = 0.0; x < 100; x += dx){
				auto val = exp(exp_mean_probability(x,mean));
				sum_st.push_back(sum);
				sum += val*dx;
			}
			
			auto k = 0u;
			for(auto x = 0.0; x < 100; x += dx){
				if(k%1000 == 0){
					cout << x << " " << log(1-sum_st[k]/sum) << " " << exp_mean_upper_probability(x,mean) << endl;
				}
				k++;
			}	
		}
		break;
	}
	
	vector <double> store;
	auto shape = 2.2, scale = 1.5, rate = 0.5;
	
	string warn;
	
	vector <double> tim_store;
	for(auto type = 0u; type < 15; type++){
		auto cl_store = clock();
		for(auto loop = 0u; loop < 10000000; loop++){
			switch(type){
			case 0: store.push_back(lognormal_sample(mean,cv,warn)); break;
			case 1: store.push_back(lognormal_probability(3,mean,cv)); break;
			case 2: store.push_back(lognormal_upper_probability(3,mean,cv)); break;
			
			case 3: store.push_back(weibull_sample(scale,shape,warn)); break;
			case 4: store.push_back(weibull_probability(3,scale,shape)); break;
			case 5: store.push_back(weibull_upper_probability(3,scale,shape)); break;
			
			case 6: store.push_back(gamma_sample(mean,cv,warn)); break;
			case 7: store.push_back(gamma_probability(3,mean,cv)); break;
			case 8: store.push_back(gamma_upper_probability(3,mean,cv)); break;
			
			case 9: store.push_back(exp_rate_sample(rate,warn)); break;
			case 10: store.push_back(exp_rate_probability(3,rate)); break;
			case 11: store.push_back(exp_rate_upper_probability(3,rate)); break;
			
			case 12: store.push_back(exp_mean_sample(mean,warn)); break;
			case 13: store.push_back(exp_mean_probability(3,mean)); break;
			case 14: store.push_back(exp_mean_upper_probability(3,mean)); break;
			}
		}
		tim_store.push_back(clock()-cl_store);
	}
		
	auto max = 0.0; for(auto val : tim_store) if(val > max) max = val;
	
	for(auto type = 0u; type < 15; type++){
		cout << "// ";
		switch(type){
		case 0: cout << "lognormal sample"; break;
		case 1: cout << "lognormal probability"; break;
		case 2: cout << "lognormal upper probability"; break;
		case 3: cout << "weibull sample"; break;
		case 4: cout << "weibull probability"; break;
		case 5: cout << "weibull upper probability"; break;
		case 6: cout << "gamma sample"; break;
		case 7: cout << "gamma probability"; break;
		case 8: cout << "gamma upper probability"; break;
		case 9: cout << "exp_rate sample"; break;
		case 10: cout << "exp_rate probability"; break;
		case 11: cout << "exp_rate upper probability"; break;
		case 12: cout << "exp_mean sample"; break;
		case 13: cout << "exp_mean probability"; break;
		case 14: cout << "exp_mean upper probability"; break;
		}
		cout << " - " << int(100.0*tim_store[type]/max) << "%" << endl;
	}
}


/// Gets M from K
vector < vector <double> > M_from_K(vector < vector <double> > &K)
{
	auto M = K;	
	
	auto N = M.size();
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			if(j != i){
				M[j][i] *= sqrt(M[i][i]*M[j][j]);
			}
		}
	}
	
	return M;
}


/// Checks the mvn jeffreys prior
void mvn_jeffreys_check()
{
	auto N = 3u;
	
	vector < vector <double> > M;

	M.resize(N);
	for(auto j = 0u; j < N; j++){
		M[j].resize(N);
		for(auto i = 0u; i < N; i++){
			if(j == i) M[j][i] = 1;
			else M[j][i] = 0;
		}
	}
	
	auto det = determinant_fast(M);
	
	auto nac = 0.0, ntr = 0.0;
	
	ofstream fout("Output/cor.txt");
	
	fout << "State";
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i <= j; i++){
			fout << "\tM" << j << i;
		}
	}
	fout << "\t" << det << endl;

	string warn;

	auto var_min = 0.5, var_max = 4.0;
	//auto var_min = 0.2, var_max = 4.0;

	//auto det_min = 0.33*N*log(var_min);
	auto det_min = N*log(var_min);

	auto smax = 1000000u;

	for(auto s = 0u;  s < smax; s++){
		if(s%1000 == 0){
			fout << s;
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i <= j; i++){
					if(j == i) fout << "\t" << M[j][i];
					else fout << "\t" << M[j][i]/sqrt(M[i][i]* M[j][j]);
				}
			}			
			fout << "\t" << det << endl;
		}
		
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i <= j; i++){
				auto d = normal_sample(0,0.1,warn);
				M[j][i] += d; M[i][j] += d;
				auto det_new = determinant_fast(M);
				ntr++;
				//if(det_new == UNSET || det_new < det_min || (i == j && M[i][i] > var_max) || (i == j && M[i][i] < var_min)){
				if(det_new == UNSET || det_new < det_min || (i == j && M[i][i] > var_max)){
					M[j][i] -= d; M[i][j] -= d;
				}
				else{
					auto al = exp(-(0.5*N+0.5)*(det_new-det));
					if(ran() < al){
						nac++;
						det = det_new;
					}
					else{
						M[j][i] -= d; M[i][j] -= d;
					}
				}
			}
		}
	}		
	cout << nac/ntr << " acc" << endl;
	
	
	vector < vector <double> > K;

	K.resize(N);
	for(auto j = 0u; j < N; j++){
		K[j].resize(N);
		for(auto i = 0u; i < N; i++){
			if(j == i) K[j][i] = 1;
			else K[j][i] = 0;
		}
	}
	
	M = M_from_K(K);
	
	det = determinant_fast(M);
	auto diag = diag_log_sum(M);
				
	ofstream fout2("Output/cor2.txt");
	
	fout2 << "State";
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i <= j; i++){
			fout2 << "\tM" << j << i;
		}
	}
	fout2 << "\t" << det << endl;

	auto si = 1.0, si2 = 1.0;
	
	auto burn = smax/10;
	
	vector < vector <double> > nacc, ntrr;
	nacc.resize(N); ntrr.resize(N);
	for(auto j = 0u; j < N; j++){
		nacc[j].resize(N,0); ntrr[j].resize(N,0);
	}
	
	for(auto s = 0u;  s < smax; s++){
		if(s%1000 == 0){
			fout2 << s;
			for(auto j = 0u; j < N; j++){
				for(auto i = 0u; i <= j; i++){
					fout2 << "\t" << K[j][i];
				}
			}			
			fout2 << "\t" << det << endl;
		}
		
		for(auto j = 0u; j < N; j++){
			for(auto i = 0u; i <= j; i++){
				double d;
				if(i == j) d = normal_sample(0,si,warn);
				else d = normal_sample(0,si2,warn);
				
				K[j][i] += d; K[i][j] += d;
				M = M_from_K(K);
				
				auto det_new = determinant_fast(M);
				auto diag_new = diag_log_sum(M);
				
				ntrr[j][i]++;
				//if(det_new == UNSET || det_new < det_min || (i == j && M[i][i] > var_max) || (i == j && M[i][i] < var_min)){
				if(det_new == UNSET || det_new < det_min || (i == j && M[i][i] > var_max)){
					K[j][i] -= d; K[i][j] -= d;
					if(s < burn){
						if(i == j) si *= 0.999;
						else si2 *= 0.999;
					}
				}
				else{
					//auto al = exp(-(0.5*N+0.5)*(det_new-det));
					auto al = exp((0.5*N-0.5)*(diag_new-diag) - (0.5*N+0.5)*(det_new-det));
					if(ran() < al){
						nacc[j][i]++;
						det = det_new;
						diag = diag_new;
						if(s < burn){
							if(i == j) si *= 1.001;
							else si2 *= 1.001;
						}
					}
					else{
						K[j][i] -= d; K[i][j] -= d;
						if(s < burn){
							if(i == j) si *= 0.999;
							else si2 *= 0.999;
						}
					}
				}
			
			}
		}

	}		
	//cout	<< nac/ntr << " acc" << endl;
	
	for(auto j = 0u; j < N; j++){
		for(auto i = 0u; i < N; i++){
			cout << nacc[j][i]/ntrr[j][i] << ",";
		}
		cout << endl;
	}
	cout << si << " " << si2 << " jj" << endl;
}
	
	
