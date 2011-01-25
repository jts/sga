//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldSearch - Find walks through a scaffold graph
//
#ifndef SCAFFOLDSEARCH_H
#define SCAFFOLDSEARCH_H

#include "GraphSearchTree.h"
#include "ScaffoldWalk.h"

//
struct ScaffoldWalkBuilder
{
    public:
        ScaffoldWalkBuilder(ScaffoldWalkVector& outWalks);
        ~ScaffoldWalkBuilder();

        // These three functions must be provided by the builder object
        // the generic graph code calls these to describe the walks through
        // the graph
        void startNewWalk(ScaffoldVertex* pStartVertex);
        void addEdge(ScaffoldEdge* pEdge);
        void finishCurrentWalk();

    private:
        ScaffoldWalkVector m_outWalks;
        ScaffoldWalk* m_pCurrWalk;
};

//
namespace ScaffoldSearch
{
    void findVariantWalks(ScaffoldVertex* pX, 
                          EdgeDir initialDir, 
                          int maxDistance,
                          size_t maxWalks, 
                          ScaffoldWalkVector& outWalks);

};

#endif
