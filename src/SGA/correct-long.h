//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// correct-long - Correct sequencing errors in long reads
//
#ifndef CORRECTLONG_H
#define CORRECTLONG_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int correctLongMain(int argc, char** argv);

// options
void parseCorrectLongOptions(int argc, char** argv);

#endif
