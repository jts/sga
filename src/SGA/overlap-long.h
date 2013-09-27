//-----------------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
// Released under the GPL
//-----------------------------------------------------
//
// overlap-long - compute overlaps for long reads
//
#ifndef OVERLAP_LONG_H
#define OVERLAP_LONG_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int overlapLongMain(int argc, char** argv);

// options
void parseOverlapLongOptions(int argc, char** argv);

#endif
