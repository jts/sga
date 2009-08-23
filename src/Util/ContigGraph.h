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

//
// Typedefs
//
typedef Bigraph ContigGraph;

void loadVertices(ContigGraph& graph, int /*kmer*/, std::string filename);
void loadEdges(ContigGraph& graph, int overlap, std::string filename);
ContigGraph* createContigGraph(int k, std::string contigsCAF, std::string adjCAF);

#endif
