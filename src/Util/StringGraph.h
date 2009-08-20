//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// String Graph - Bidirectional graph of sequence reads
// and their overlaps
//
#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H

#include "Bigraph.h"
#include "Contig.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>

typedef SeqItem StringNode;
//
// Typedefs
//
typedef Edge<Sequence> StringEdge;
typedef Vertex<StringNode, StringEdge> StringVertex;
typedef Bigraph<StringVertex> StringGraph;

StringGraph* createStringGraph(std::string readFile, std::string overlapFile);

#endif
