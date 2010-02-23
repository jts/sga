//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Stats - Common statistics functions used 
// throughout the program
//

#ifndef STATS_H
#define STATS_H

namespace Stats
{
//
// Probability
//
double logPoisson(unsigned int k, double m);
double logFactorial(unsigned int k);

};
#endif
