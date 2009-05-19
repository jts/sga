#include <iostream>
#include <math.h>
#include "Util.h"

//
// Globals
//

//
// Sequence operations
//

//
// Reverse complement a sequence
//
Sequence reverseComplement(Sequence seq)
{
	std::string out(seq.length(), 'A');
	size_t last_pos = seq.length() - 1;
	for(int i = last_pos; i >= 0; --i)
	{
		out[last_pos - i] = complement(seq[i]);
	}
	return out;
}

//
// Complement a base
//
char complement(char base)
{
	switch(base)
	{
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		default:
			assert(false && "Unknown base!");
	}
}

//
// Probability functions
//

//
// Poission 
//
double poisson(int k, double m)
{

	double f_k = factorial(k);
	double p = pow(m, k) * exp(-m);
	//std::cout << "k: " << k << " f: " << f_k << " m: " << m << " p: " << p << std::endl;
	return p / f_k;
}

//
// Factorial
//
int factorial(int k)
{
	int result = 1;
	while(k > 0)
		result *= k--;
	return result;
}

//
// Poisson in log space
//
double log_poisson(int k, double m)
{
	double f_k = log_factorial(k);
	double p = (double)k * log(m) - m - f_k;
	//std::cout << "k: " << k << " f: " << f_k << " m: " << m << " p: " << p << std::endl;
	return p;
}

//
// Factorial in log space
//
double log_factorial(int k)
{
	double result = 0;
	while(k > 0)
		result += log(k--); //slow
	return result;
}


