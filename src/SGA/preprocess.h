//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// preprocess - prepare data files for assembly
//
#ifndef PREPROCESS_H
#define PREPROCESS_H
#include <getopt.h>
#include "config.h"
#include "Quality.h"

// functions
int preprocessMain(int argc, char** argv);
void parsePreprocessOptions(int argc, char** argv);
bool processRead(SeqRecord& record);
bool samplePass();
char randomBase();
void softClip(int qualTrim, std::string& seq, std::string& qual);

#endif
