//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// rmdup - Remove duplicate reads from the data set
//
#ifndef RMDUP_H
#define RMDUP_H
#include <getopt.h>
#include "config.h"

// functions
int rmdupMain(int argc, char** argv);
void parseRmdupOptions(int argc, char** argv);
void rmdup();
void parseDupHits(const StringVector& hitsFilenames);

#endif
