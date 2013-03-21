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
