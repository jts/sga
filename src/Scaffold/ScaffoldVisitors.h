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
                                             int& dist, double& sd, EdgeDir& dir_yz, 
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
        size_t m_sumUnique;
        size_t m_sumRepeat;
};

// Walk the graph validating that the links present are valid
class ScaffoldLinkValidator
{
    public:
        
        ScaffoldLinkValidator(int maxOverlap, double threshold);
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

// Detect and remove polymorphmic vertices in the scaffold
class ScaffoldPolymorphismVisitor
{
    public:
        ScaffoldPolymorphismVisitor();
        
        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/);
};

// Walk through the graph breaking up components whenever a vertex has more
// than one edge in a particular direction
class ScaffoldMultiEdgeRemoveVisitor
{
    public:
        
        ScaffoldMultiEdgeRemoveVisitor() {}
        void previsit(ScaffoldGraph* /*pGraph*/) {}
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/) {}

};

// Output scaffolds
class ScaffoldWriterVisitor
{
    public:
        ScaffoldWriterVisitor(const std::string& filename);
        ~ScaffoldWriterVisitor();

        void previsit(ScaffoldGraph* /*pGraph*/);
        bool visit(ScaffoldGraph* pGraph, ScaffoldVertex* pVertex);
        void postvisit(ScaffoldGraph* /*pGraph*/) {}

    private:
        std::ostream* m_pWriter;
};

#endif
