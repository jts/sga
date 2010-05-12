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

// Add overlaps to pX inferred from the edges of pY to outMap
void addOverlapsToSet(const Vertex* pX, const EdgeDesc& edXY, const Overlap& ovrXY, 
                      double maxER, int minLength, bool bothDirs, EdgeDescOverlapMap& outMap);

// Walk the graph to collect the entire set of overlaps for pVertex via its neighbors
void constructCompleteOverlapMap(const Vertex* pVertex, double maxER, int minLength, 
                                 bool bothDirs, EdgeDescOverlapMap& outMap);

// Partition the complete overlap set of pVertex into irreducible and transitive edge sets
void constructPartitionedOverlapMap(const Vertex* pVertex, double maxER, int minLength, 
                                    double reduceER, int reduceLength,
                                    EdgeDescOverlapMap& irreducibleMap, 
                                    EdgeDescOverlapMap& transitiveMap);

// Move transitive overlaps in irreducibleMap to the transitiveMap
void partitionOverlapMap(double maxER, int minLength, 
                         EdgeDescOverlapMap& irreducibleMap, 
                         EdgeDescOverlapMap& transitiveMap);

// Diff the overlaps indicated by inMap with the edge list of a vertex
// The edges that are in inMap by not pVertex are placed into missingMap
// The edges that are in pVertex but not inMap are added to extraMap
void diffVertexOverlapMap(Vertex* pVertex, const EdgeDescOverlapMap& inMap, 
                          EdgeDescOverlapMap& missingMap,
                          EdgeDescOverlapMap& extraMap);

// Filter out overlaps from an overlapMap 
void filterOverlapMap(const Vertex* pVertex, double maxER, int minLength, EdgeDescOverlapMap& overlapMap);

// Calculate the error rate between the two vertices
double calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY);

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

//
// Construct an extended multioverlap for a vertex
//
MultiOverlap makeExtendedMultiOverlap(const Vertex* pVertex);

//
// Construct SeqTries from the extended overlap set
//
void makeExtendedSeqTries(const Vertex* pVertex, double p_error, SeqTrie* pLeftTrie, SeqTrie* pRightTrie);

// Simple getters for std::transform
EdgeDesc getEdgeDescFromEdge(Edge* pEdge);
EdgeDesc getEdgeDescFromPair(const EdgeDescOverlapPair& pair);

void printOverlapMap(const EdgeDescOverlapMap& overlapMap);

};

#endif
