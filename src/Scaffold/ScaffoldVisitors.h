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

//
namespace ScaffoldAlgorithms
{
    bool areEdgesAmbiguous(ScaffoldEdge* pXY, ScaffoldEdge* pXZ);
    void inferScaffoldEdgeYZ(ScaffoldEdge* pXY, ScaffoldEdge* pXZ,
                                             int& dist, EdgeDir& dir_yz, 
                                             EdgeDir& dir_zy, EdgeComp& comp);

};

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

// Classify vertices as repetitive based on their set of edges
class ScaffoldEdgeSetClassificationVisitor
{
    public:
        
        ScaffoldEdgeSetClassificationVisitor(int maxOverlap, double threshold);
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:

        int m_maxOverlap;
        double m_threshold;
        size_t m_numUnique;
        size_t m_numRepeat;
};

// Walk through the scaffold graph making chains of vertices
class ScaffoldChainVisitor
{
    public:
        
        ScaffoldChainVisitor(int maxOverlap);
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);

    private:
        int m_maxOverlap;

};
#endif
