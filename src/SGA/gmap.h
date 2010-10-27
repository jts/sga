//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// gmap - Map sequences to the vertices of a graph
//
#ifndef GMAP_H
#define GMAP_H
#include <getopt.h>
#include "config.h"

// functions
int gmapMain(int argc, char** argv);
void parseGmapOptions(int argc, char** argv);
void gmap();
void parseGmapHits(const StringVector& hitsFilenames);

#endif
