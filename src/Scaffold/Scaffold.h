#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include "Util.h"
#include "Contig.h"
#include "Bigraph.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>

struct ScaffoldData
{
    ScaffoldData(int d, double sd, int np) : estDist(d), stdDev(sd), numPairs(np) {}
    int estDist;
    double stdDev;
    unsigned int numPairs;

    friend std::ostream& operator<<(std::ostream& out, const ScaffoldData& sd)
    {
        out << sd.estDist;
        return out;
    }
};

struct SLink
{
    ContigID linkedID;
    int dist;
    unsigned int numPairs;
    double stdDev;
    bool isRC;
    
    friend std::istream& operator>>(std::istream& in, SLink& sl)
    {
        std::string line;
        in >> line;
        if(line.size () != 0)
        {
            StringVec fields = split(line, ',');
            assert(fields.size() == 5);
            sl.linkedID = fields[0];
            sl.dist = atoi(fields[1].c_str());
            sl.numPairs = atoi(fields[2].c_str());
            std::stringstream fparser(fields[3]);
            fparser >> sl.stdDev;
            sl.isRC = atoi(fields[4].c_str());
        }
        return in;
    }

    friend std::ostream& operator<<(std::ostream& out, const SLink& sl)
    {
        out << sl.linkedID << " " << sl.dist << " " << sl.numPairs << " " << sl.stdDev << " " << sl.isRC;
        return out;
    }        
};


//
// Typedefs
//
typedef Edge<ScaffoldData> ScaffoldEdge;
typedef Vertex<Contig, ScaffoldEdge> ScaffoldVertex;
typedef Bigraph<ScaffoldVertex> ScaffoldGraph;


//
// Structs
//
struct LinearScaffoldLink
{
    LinearScaffoldLink(ScaffoldEdge e, Range r) : edge(e), range(r) {}
    ScaffoldEdge edge;
    Range range;

    friend std::ostream& operator<<(std::ostream& out, const LinearScaffoldLink& lsl)
    {
        out << lsl.edge.getEnd() << " " << lsl.range;
        return out;
    }

    static bool sortStarts(const LinearScaffoldLink& lsl1, const LinearScaffoldLink& lsl2)
    {
        return lsl1.range.start < lsl2.range.start;
    }
};

//
// Typedefs
//
typedef std::map<ContigID, ScaffoldData> SDMap;
typedef std::vector<Range> RangeVec;
typedef std::vector<LinearScaffoldLink> LSLVec;

//
// Functions
//
void buildGraph(SDMap& sdMap);
bool areEdgesConsistent(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex, EdgeDir dir);

// Vertex visit functions
bool makeTransitive(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
bool cutInconsistent(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
bool cutAmbigious(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);

Contig& getContig(SDMap& sdMap, ContigID cid);
Range convertEdgeToRange(const ScaffoldGraph* sg, const ScaffoldEdge& e);

// Graph building
void addVertexToScaffoldGraph(ScaffoldGraph& graph, SDMap& sdMap, ContigID id);
void addEdgeToScaffoldGraph(ScaffoldGraph& graph, ContigID id1, ContigID id2, EdgeDir dir, EdgeComp comp, int dist);

// Parsing
void parseLinks(std::string filename, ScaffoldGraph& graph);
void parseOptions(int argc, char** argv);

// Writing
void writeScaffold(ostream& out, int idNum, const ScaffoldGraph::Path& path);
void writeScaffoldNode(ostream& out, VertexID id, int dist, bool orientation);

//
// Getopt
//
#define PROGRAM "Scaffold"
static const char *VERSION_MESSAGE =
PROGRAM "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION] ... CONTIGS DISTANCES\n"
"Build scaffolds from the CONTIGS based on the inferred DISTANCES.\n"
"\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n";


namespace opt
{
    static unsigned int verbose;
}

static const char* shortopts = "o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


#endif
