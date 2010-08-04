//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVisitors - Functors that perform
// some operation on a ScaffoldVertex/Graph
//
#ifndef SCAFFOLDVISITORS_H
#define SCAFFOLDVISITORS_H

#include "ScaffoldGraph.h"

class ScaffoldStatsVisitor
{
    public:
        
        void previsit(ScaffoldGraph* pGraph);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* pGraph);

    private:
        size_t m_numVertices;
};

#endif
