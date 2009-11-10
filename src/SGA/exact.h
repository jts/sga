//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// exact - Test algorithms for assembling exact data
//
#ifndef EXACT_H
#define EXACT_H
#include <getopt.h>
#include "config.h"
#include "SGACommon.h"

// functions
int exactMain(int argc, char** argv);
void parseExactOptions(int argc, char** argv);
void exact();

#endif
