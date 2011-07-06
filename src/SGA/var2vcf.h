//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// var2vcf - convert variants from sga graph-diff
// or sga assembly to a vcf file of differences
// with respect to a reference
//
#ifndef VAR2VCF_H
#define VAR2VCF_H
#include <getopt.h>
#include "config.h"
#include "BWT.h"
#include "Match.h"
#include "BWTAlgorithms.h"
#include "OverlapAlgorithm.h"

// functions

//
int var2vcfMain(int argc, char** argv);

// options
void parseVar2VCFOptions(int argc, char** argv);

#endif
