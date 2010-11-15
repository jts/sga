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
#include <deque>

// A walk on the string graph is given by the starting vertex
// then a vector of edges used in the walk
enum SGWalkType
{
    SGWT_START_TO_END,
    SGWT_EXTENSION
};

class SGWalk
{
    public:
        
        SGWalk(Vertex* pStartVertex, bool bIndexWalk = false);
        SGWalk(const SGWalk& other);

        ~SGWalk();

        SGWalk& operator=(const SGWalk& other);


        void addEdge(Edge* pEdge);
        void popLast();

        VertexPtrVec getVertices() const;
        Vertex* getStartVertex() const;
        Edge* getLastEdge() const;
        Edge* getEdge(size_t idx) const;
        size_t getNumVertices() const;
        size_t getNumEdges() const;

        // Returns true if the walk contains the specified vertex
        // If the walk is not indexed, this will assert
        bool containsVertex(const VertexID& id) const;

        // Truncate the walk after the first instance of id
        void truncate(const VertexID& id);

        //
        std::string getString(SGWalkType type) const;

        // distance calculations
        int getExtensionDistance() const;
        int getEndToEndDistance() const;
        int getStartToEndDistance() const;
        int getEndToStartDistance() const;

        void print() const;

    private:
        
        Vertex* m_pStartVertex;
        EdgePtrVec m_edges;
        
        typedef std::set<VertexID> WalkIndex;
        WalkIndex* m_pWalkIndex;

        // The distance from the end of pStart to the last vertex in the walk
        // This equals the length of the extension string
        // x -----------
        // y     ------------
        // z              ----------
        // distance     ************
        int m_extensionDistance;
};
typedef std::vector<SGWalk> SGWalkVector;
typedef std::deque<SGWalk> WalkQueue;

// String Graph searching algorithms
namespace SGSearch
{
    //
    void findWalks(Vertex* pX, Vertex* pY, EdgeDir initialDir,
                   int maxDistance, size_t maxQueue, SGWalkVector& outWalks);

    void findCollapsedWalks(Vertex* pX, EdgeDir initialDir, 
                            int maxDistance, size_t maxQueue, 
                            SGWalkVector& outWalks);

    // Count the number of vertices that span the sequence junction
    // described by edge XY. Returns -1 if the search was not completed
    int countSpanningCoverage(Edge* pXY, size_t maxQueue);

    //
    void initializeWalkQueue(Vertex* pX, EdgeDir initialDir, bool bIndexWalks, WalkQueue& queue);
    bool extendWalk(const Vertex* pX, EdgeDir dir, SGWalk& currWalk, WalkQueue& queue);
};

#endif
