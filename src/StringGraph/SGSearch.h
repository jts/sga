//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGSearch - Algorithms and data structures
// for searching a string graph
//
#ifndef SGSEARCH_H
#define SGSEARCH_H

#include "Bigraph.h"
#include "SGWalk.h"
#include <deque>

// String Graph searching algorithms
namespace SGSearch
{
    //
    void findWalks(Vertex* pX, Vertex* pY, EdgeDir initialDir,
                   int maxDistance, size_t maxNodes, SGWalkVector& outWalks);

    void findVariantWalks(Vertex* pX, 
                          EdgeDir initialDir, 
                          int maxDistance,
                          size_t maxWalks, 
                          SGWalkVector& outWalks);

    void findCollapsedWalks(Vertex* pX, EdgeDir initialDir, 
                            int maxDistance, size_t maxNodes,
                            SGWalkVector& outWalks);

    void convertEdgeVectorsToSGWalk(const std::vector<EdgePtrVec>& edgeWalks, bool bIndexWalks, SGWalkVector& outWalks);

    // Count the number of vertices that span the sequence junction
    // described by edge XY. Returns -1 if the search was not completed
    int countSpanningCoverage(Edge* pXY, size_t maxQueue);

    // Returns true if all the endpoints of the edges in epv are in vertexSet
    bool checkEndpointsInSet(EdgePtrVec& epv, std::set<Vertex*>& vertexSet);
};

#endif
