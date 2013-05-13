//-----------------------------------------------
// Copyright 2009-2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAStats - Common statistics functions used 
// throughout the program
//
#include "SGAStats.h"
#include <math.h>
#include <iostream>
#include <assert.h>

//
// Poisson in log space
//
double SGAStats::logPoisson(unsigned int k, double m)
{
    double f_k = logFactorial(k);
    double p = (double)k * log(m) - m - f_k;
    //std::cout << "k: " << k << " f: " << f_k << " m: " << m << " p: " << p << std::endl;
    return p;
}

//
// Factorial in log space
//
double SGAStats::logFactorial(unsigned int k)
{
    double result = 0;
    while(k > 0)
        result += log(k--); //slow
    return result;
}

// 
// Log binomial pmf 
//
double SGAStats::logBinomial(unsigned int k, unsigned int n, double p)
{
    assert(k <= n);
    double b_coeff = logFactorial(n) - logFactorial(k) - logFactorial(n - k);
    double lp = k * log(p) + (n - k) * log(1 - p);
    return b_coeff + lp; 
}

//
//
//
double SGAStats::logIntegerBeta(double x, unsigned int a, unsigned int b)
{
    assert(x >= 0.0f && x <= 1.0f);
    assert(a > 0 && b > 0);

    double log_beta_f = logFactorial(a - 1) + logFactorial(b - 1) - logFactorial(a + b - 1);
    double log_p = (a - 1) * log(x) + (b - 1) * log(1 - x) - log_beta_f;
    return log_p;
}

