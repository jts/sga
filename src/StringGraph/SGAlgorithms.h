//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAlgorithms - Algorithms for manipulating
// the string graph
//
#ifndef SGALGORITHMS_H
#define SGALGORITHMS_H

#include "Bigraph.h"
#include "SGUtil.h"
#include <queue>

namespace SGAlgorithms
{

//
// Helper data structures
//

typedef std::pair<EdgeDesc, Overlap> EdgeDescOverlapPair;

// Comparator
struct EDOPairCompare
{
    bool operator()(const EdgeDescOverlapPair& edpXY, const EdgeDescOverlapPair& edpXZ)
    {
        return edpXY.second.match.coord[0].length() < edpXZ.second.match.coord[0].length();
    }
};

// typedefs
typedef std::priority_queue<EdgeDescOverlapPair, 
                            std::vector<EdgeDescOverlapPair>,
                            EDOPairCompare> EDOPairQueue;

typedef std::map<EdgeDesc, Overlap> EdgeDescOverlapMap;
typedef std::set<EdgeDesc> EdgeDescSet;

// Find new edges for pVertex that are required if pDeleteEdge is removed from the graph
void remodelVertexForExcision(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge);

// Create the edges described by the overlap.
Edge* createEdgesFromOverlap(StringGraph* pGraph, const Overlap& o, bool allowContained);

// Calculate the error rate between the two vertices
double calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY);

void updateContainFlags(StringGraph* pGraph, Vertex* pVertex, EdgeDescOverlapMap& containMap);
void updateContainFlags(StringGraph* pGraph, Vertex* pVertex, const EdgeDesc& ed, const Overlap& ovr);

//
// Overlap inference algorithms
//
// Infer an overlap from two edges
// The input overlaps are between X->Y Y->Z
// and the returned overlap is X->Z
Overlap inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);
EdgeDesc overlapToEdgeDesc(Vertex* pY, const Overlap& ovrXY);

// Returns true if XZ has a non-zero length overlap
bool hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);

// Construct an extended multioverlap for a vertex
MultiOverlap makeExtendedMultiOverlap(const StringGraph* pGraph, const Vertex* pVertex);

// Construct SeqTries from the extended overlap set
void makeExtendedSeqTries(const StringGraph* pGraph, const Vertex* pVertex, 
                          double p_error, SeqTrie* pLeftTrie, SeqTrie* pRightTrie);

// Simple getters for std::transform
EdgeDesc getEdgeDescFromEdge(Edge* pEdge);
EdgeDesc getEdgeDescFromPair(const EdgeDescOverlapPair& pair);

void printOverlapMap(const EdgeDescOverlapMap& overlapMap);

};

#endif
