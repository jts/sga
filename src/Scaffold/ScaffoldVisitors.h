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

// Write summary statistics of the scaffold graph to stdout
class ScaffoldStatsVisitor
{
    public:
        
        void previsit(ScaffoldGraph* pGraph);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* pGraph);

    private:
        size_t m_numVertices;
        size_t m_numEdges;
};

// Classify vertices as unique or repetitive based on their
// a-statistic values
class ScaffoldAStatisticVisitor
{
    public:
        
        ScaffoldAStatisticVisitor(double uniqueThreshold, double repeatThreshold);
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        double m_uniqueThreshold;
        double m_repeatThreshold;

        size_t m_numUnique;
        size_t m_numRepeat;
};

#endif
