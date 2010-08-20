///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringGraphGenerator - Iteratively construct 
// a local String Graph using the FM-index
//
#ifndef STRINGGRAPHGENERATOR_H
#define STRINGGRAPHGENERATOR_H

#include "OverlapAlgorithm.h"
#include "SGUtil.h"
#include "GraphCommon.h"
#include <queue>

// The GraphFrontier is a node that is on the edge
// of the graph - it can be used to search the FM-index
// to expand the graph outwards.
struct GraphFrontier
{
    Vertex* pVertex;
    EdgeDir dir;
    int distance; // The distance from the starting node
};

typedef std::queue<GraphFrontier> FrontierQueue;

class StringGraphGenerator
{
    public:
        StringGraphGenerator(const OverlapAlgorithm* pOverlapper,
                             const SeqRecord& startRead, 
                             const SeqRecord& endRead, 
                             int minOverlap,
                             EdgeDir startDir);

        ~StringGraphGenerator();


    private:

        //
        void buildGraph(FrontierQueue& queue, int maxDistance);
        void updateGraphAndQueue(GraphFrontier& currNode, FrontierQueue& queue, OverlapBlockList& blockList, int maxDistance);
        
        Vertex* addTerminalVertex(const SeqRecord& record);

        //
        std::string overlapBlockToCanonicalID(OverlapBlock& block);
        
        // Data
        const OverlapAlgorithm* m_pOverlapper;
        int m_minOverlap;
        StringGraph* m_pGraph;

        Vertex* m_pStartVertex;
        Vertex* m_pEndVertex;

        static const GraphColor UNEXPLORED_COLOR = GC_WHITE;
        static const GraphColor EXPLORED_COLOR = GC_BLACK;
};

#endif
