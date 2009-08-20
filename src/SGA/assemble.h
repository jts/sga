//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// assemble - Assemble reads into contigs
//
#ifndef ASSEMBLE_H
#define ASSEMBLE_H
#include <getopt.h>
#include "config.h"

// functions
int assembleMain(int argc, char** argv);
void parseAssembleOptions(int argc, char** argv);
void assemble();

#endif
