#include "StatsCommon.h"
#include <math.h>
#include <iostream>

//
// Probability functions
//

//
// Poisson in log space
//
double log_poisson(unsigned int k, double m)
{
	double f_k = log_factorial(k);
	double p = (double)k * log(m) - m - f_k;
	//std::cout << "k: " << k << " f: " << f_k << " m: " << m << " p: " << p << std::endl;
	return p;
}

//
// Factorial in log space
//
double log_factorial(unsigned int k)
{
	double result = 0;
	while(k > 0)
		result += log(k--); //slow
	return result;
}

