//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// String Graph - Bidirectional graph of sequence reads
// and their overlaps
//
#ifndef STRINGGRAPH_H
#define STRINGGRAPH_H

#include "Match.h"
#include "Bigraph.h"
#include <boost/pool/object_pool.hpp>
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>

#define SE_CAST(x) static_cast<StringEdge*>((x))
#define CSE_CAST(x)  static_cast<const StringEdge*>((x))
#define SV_CAST(x) static_cast<StringVertex*>((x))
#define CSV_CAST(x)  static_cast<const StringVertex*>((x))


typedef Bigraph StringGraph;
class StringVertex;
class StringEdge;

#endif

