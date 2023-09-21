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

/// Checks that a distribution is being generated correctly
void test_distribution()
{
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
				auto val = normal_sample(meansamp,sdsamp);
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
				auto val = gamma_sample(meansamp,cvsamp);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 100; x+=dx){
				auto val = exp(gamma_probability(x,meansamp,cvsamp));
				//cout << x << " " << val << "val\n";
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
				auto val = lognormal_sample(meansamp,cvsamp);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 100; x+=dx){
				auto val = exp(lognormal_probability(x,meansamp,cvsamp));
				//cout << x << " " << val << "val\n";
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
				auto val = weibull_sample(shapesamp,scalesamp);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 100; x+=dx){
				auto val = exp(weibull_probability(x,shapesamp,scalesamp));
				//cout << x << " " << val << "val\n";
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
				auto val = poisson_sample(lamsamp);
				sum += val; sum2 += val*val;
			}

			for(auto k = 0u; k < 100; k++){
				auto val = exp(poisson_probability(k,lamsamp));
				//cout << x << " " << val << "val\n";
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
				auto val = bernoulli_sample(zsamp);
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
				auto val = beta_sample(alphasamp,betasamp);
				sum += val; sum2 += val*val;
			}

			for(auto x = dx; x < 1-dx; x+=dx){
				auto val = exp(beta_probability(x,alphasamp,betasamp));
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
