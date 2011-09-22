//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// metagenome - Assemble metagenomics data
//
#ifndef METAGENOME_H
#define METAGENOME_H
#include <getopt.h>
#include "config.h"

// functions

//
int metagenomeMain(int argc, char** argv);

// options
void parseMetagenomeOptions(int argc, char** argv);

#endif
