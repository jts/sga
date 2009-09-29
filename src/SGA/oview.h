//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// oview - view overlaps between reads
//
#ifndef OVIEW_H
#define OVIEW_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"

// functions

//
int oviewMain(int argc, char** argv);

// options
void parseOviewOptions(int argc, char** argv);

#endif
