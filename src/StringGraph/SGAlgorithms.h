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
// Overlap discovery algorithms
//

// Structure used for iterative exploration of the graph
// The exploration starts at some vertex X, each element
// holds a possible overlap between X and some other vertex Y
struct ExploreElement
{
    ExploreElement(const EdgeDesc& e, const Overlap& o) : ed(e), ovr(o) {}
    EdgeDesc ed;
    Overlap ovr;
};

// Comparison operator used to compare ExploreElements
// by the length of the overlap on vertex X
struct CompareExploreElemOverlapLength
{
    bool operator()(const ExploreElement& elemXY, const ExploreElement& elemXZ)
    {
        return elemXY.ovr.match.coord[0].length() < elemXZ.ovr.match.coord[0].length();
    }
};

typedef std::queue<ExploreElement> ExploreQueue;
typedef std::priority_queue<ExploreElement, 
                            std::vector<ExploreElement>, 
                            CompareExploreElemOverlapLength> ExplorePriorityQueue;

//
typedef std::map<EdgeDesc, Overlap> EdgeDescOverlapMap;
typedef std::set<EdgeDesc> EdgeDescSet;

// Remodel the edges of pVertex by finding any new irreducible edges
// that may need to be added if the vertex at the end of pEdge is removed
// from the graph
void remodelVertexAfterRemoval(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge);

// Find new edges for pVertex that are required if pDeleteEdge is removed from the graph
void remodelVertexForExcision(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge);

// Add the neighbors of the endpoint of edXY to the explore queue if they overlap pX
// and are not present in pSeenEdges. Edges are added to pSeenEdges as they are found
// so the contents of pSeenEdges will be modified by this function.
void enqueueEdges(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, 
                  ExplorePriorityQueue& outQueue, EdgeDescSet* pSeenEdges);

// Add overlaps to pX inferred from the edges of pY to outMap
void addOverlapsToSet(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, EdgeDescOverlapMap& outMap);

// Discover the complete set of overlaps for pVertex
void findOverlapMap(const Vertex* pVertex, double maxER, int minLength, EdgeDescOverlapMap& outMap);

// recursive function to discover overlaps of pX from pY
void _discoverOverlaps(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, 
                       double maxER, int minLength, EdgeDescOverlapMap& outMap);

// Calculate the error rate between the two vertices
double calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY);

//
// Overlap inference algorithms
//
// Infer an overlap from two edges
// The input overlaps are between X->Y Y->Z
// and the returned overlap is X->Z
Overlap inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);
EdgeDesc inferTransitiveEdgeDesc(const EdgeDesc& edXY, const EdgeDesc& edYZ);

// Returns true if XZ has a non-zero length overlap
bool hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ);

//
// Construct an extended multioverlap for a vertex
//
MultiOverlap makeExtendedMultiOverlap(const Vertex* pVertex);

//
// Construct SeqTries from the extended overlap set
//
void makeExtendedSeqTries(const Vertex* pVertex, double p_error, SeqTrie* pLeftTrie, SeqTrie* pRightTrie);

};

#endif
