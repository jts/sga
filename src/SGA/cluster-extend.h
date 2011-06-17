//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// cluster-extend - extend previously-built clusters
//
#ifndef CLUSTER_EXTEND_H
#define CLUSTER_EXTEND_H
#include <getopt.h>
#include "config.h"

// functions
int clusterExtendMain(int argc, char** argv);
void parseClusterExtendOptions(int argc, char** argv);
void clusterExtend();

#endif
