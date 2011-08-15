//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// convert-beetl - Utility to convert a beetl-constructed
// bwt index into SGA's format
//
#ifndef BEETL_CONVERT_H
#define BEETL_CONVERT_H
#include <getopt.h>
#include "config.h"

int convertBeetlMain(int argc, char** argv);
void parseConvertBeetlOptions(int argc, char** argv);

#endif
