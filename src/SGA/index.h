//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// index - Read in sequences and index them in a suffix array
//
#ifndef INDEX_H
#define INDEX_H
#include <getopt.h>
#include "config.h"
int indexMain(int argc, char** argv);
void buildIndex(std::string filename);
void parseIndexOptions(int argc, char** argv);

#endif
