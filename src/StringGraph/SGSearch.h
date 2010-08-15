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

// A walk on the string graph is given by the starting vertex
// then a vector of edges used in the walk
class SGWalk
{
    public:
        SGWalk(const Vertex* pStartVertex);
        void addEdge(Edge* pEdge);
        Edge* getLastEdge() const;

        void print() const;

    private:
        
        const Vertex* m_pStartVertex;
        EdgePtrVec m_edges;
};
typedef std::vector<SGWalk> SGWalkVector;

// String Graph searching algorithms
namespace SGSearch
{
    void findWalks(const Vertex* pX, const Vertex* pY, EdgeDir initialDir,
                   int maxDistance, SGWalkVector& outWalks);

};

#endif
