//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// hapgen - Generate candidate haplotypes from
// an assembly graph
//
#ifndef HAPGEN_H
#define HAPGEN_H
#include <getopt.h>
#include "config.h"

// functions

//
int hapgenMain(int argc, char** argv);

// options
void parseHapgenOptions(int argc, char** argv);

#endif
