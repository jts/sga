//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVisitors - Functors that perform
// some operation on a ScaffoldVertex/Graph
//
#include "ScaffoldVisitors.h"
#include <limits.h>

//
// Returns true if the distance estimates for the two edges
// are close enough that the ordering cannot be resolved. Assumes
// that pXY is closer than pXZ
bool ScaffoldAlgorithms::areEdgesAmbiguous(ScaffoldEdge* pXY, ScaffoldEdge* pXZ)
{
    const int AMBIGUOUS_TOLERANCE = 2;
    assert(pXY->getDistance() <= pXZ->getDistance());
    int translated = pXZ->getDistance() - AMBIGUOUS_TOLERANCE*pXZ->getStdDev();
    if(translated < pXY->getDistance())
        return true;
    else
        return false;
}

//
// Returns the edge properties inferred from the edges X->Y and X->Z
void ScaffoldAlgorithms::inferScaffoldEdgeYZ(ScaffoldEdge* pXY, ScaffoldEdge* pXZ,
                                             int& dist, EdgeDir& dir_yz, 
                                             EdgeDir& dir_zy, EdgeComp& comp)
{
    comp = (pXY->getComp() == pXZ->getComp() ? EC_SAME : EC_REVERSE);
    dir_yz = !pXY->getTwin()->getDir(); // Opposite direction of X <- Y
    dir_zy = pXY->getTwin()->getDir();  // Same direction as X <- Z since Z is after Y
    dist = pXZ->getDistance() - (pXY->getDistance() + pXY->getEnd()->getSeqLen());
}

// ScaffoldStatsVisitor
void ScaffoldStatsVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_numVertices = 0;
    m_numEdges = 0;
}

bool ScaffoldStatsVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    ++m_numVertices;
    m_numEdges += pVertex->getNumEdges();

    return false;
}

void ScaffoldStatsVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    printf("Scaffold Stats -- Num vertices: %zu Num edges: %zu\n", m_numVertices, m_numEdges);
}

// ScaffoldAStatistic
ScaffoldAStatisticVisitor::ScaffoldAStatisticVisitor(double uniqueThreshold, 
                                                     double repeatThreshold) : m_uniqueThreshold(uniqueThreshold),
                                                                               m_repeatThreshold(repeatThreshold)
{
                                                                               
}

//
void ScaffoldAStatisticVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_numUnique = 0;
    m_numRepeat = 0;
}

//
bool ScaffoldAStatisticVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never re-classify repeats
    if(pVertex->getClassification() != SVC_REPEAT)
    {
        if(pVertex->getAStatistic() > m_uniqueThreshold)
        {
            pVertex->setClassification(SVC_UNIQUE);
            ++m_numUnique;
        }

        if(pVertex->getAStatistic() < m_repeatThreshold)
        {
            pVertex->setClassification(SVC_REPEAT);
            ++m_numRepeat;
        }
    }
    return false;   
}

//
void ScaffoldAStatisticVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    std::cerr << "A-statistic: classified " << m_numUnique << " components as unique\n";
    std::cerr << "A-statistic: classified " << m_numRepeat << " components as repeat\n"; 
}

//
ScaffoldEdgeSetClassificationVisitor::ScaffoldEdgeSetClassificationVisitor(int maxOverlap, 
                                                                           double threshold) : m_maxOverlap(maxOverlap),
                                                                                               m_threshold(threshold)
{

}

//
void ScaffoldEdgeSetClassificationVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
    m_numUnique = 0;
    m_numRepeat = 0;
}

//
bool ScaffoldEdgeSetClassificationVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never re-classify repeats
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    assert(false && "Not implemented");
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 0)
        {
            // Process the edges in order of length
            std::sort(edgeVec.begin(), edgeVec.end(), ScaffoldEdgePtrDistanceCompare);
        }
    }
    return false;   
}

//
void ScaffoldEdgeSetClassificationVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    std::cerr << "Edge set overlap: classified " << m_numRepeat << " components as repeat\n"; 
}

//
ScaffoldChainVisitor::ScaffoldChainVisitor(int maxOverlap) : m_maxOverlap(maxOverlap)
{

}

//
void ScaffoldChainVisitor::previsit(ScaffoldGraph* /*pGraph*/)
{
}

// 
bool ScaffoldChainVisitor::visit(ScaffoldGraph* /*pGraph*/, ScaffoldVertex* pVertex)
{
    // Never try to make chains from a repeat
    if(pVertex->getClassification() == SVC_REPEAT)
        return false;

    bool changed_graph = false;
    for(size_t idx = 0; idx < ED_COUNT; idx++)
    {
        EdgeDir dir = EDGE_DIRECTIONS[idx];
        ScaffoldEdgePtrVector edgeVec = pVertex->getEdges(dir);
        if(edgeVec.size() > 1)
        {
            // Try to resolve this vertex edge set
            // by creating edges between the vertices that it
            // is linked to
            
            // Process the edges in order of length
            std::sort(edgeVec.begin(), edgeVec.end(), ScaffoldEdgePtrDistanceCompare);
            ScaffoldEdge* pXY = edgeVec[0];
            ScaffoldEdge* pXZ = edgeVec[1];
            assert(pXY->getDistance() <= pXZ->getDistance());
            ScaffoldVertex* pY = pXY->getEnd();
            ScaffoldVertex* pZ = pXZ->getEnd();

            if(pY->isRepeat() || pZ->isRepeat())
            {
                std::cerr << "Warning, skipping repeat\n";
                continue;
            }

            // Check if the ordering of the contigs is ambiguous
            bool isAmbiguous = ScaffoldAlgorithms::areEdgesAmbiguous(pXY, pXZ);
            if(isAmbiguous)
            {
                std::cerr << "Warning edges " << *pXY << " and " << *pXZ << 
                              " are ambiguously ordered\n";
            }

            // Infer the edge data for the edge Y->Z
            int dist;
            EdgeDir dir_yz;
            EdgeDir dir_zy; 
            EdgeComp comp;
            ScaffoldAlgorithms::inferScaffoldEdgeYZ(pXY, pXZ, dist, dir_yz, dir_zy, comp);

            // Sanity checks
            bool isConsistent = (dist > -1*m_maxOverlap);
            if(!isConsistent)
            {
                std::cout << "Inferred edge is not consistent with max overlap: " << dist << "\n";
                std::cout << "Input: " << *pXY << " " << *pXZ << "\n";
            }

            // Check if an edge between Y->Z already exists
            ScaffoldEdge* pCheckEdge = pY->findEdgeTo(pZ->getID(), dir_yz, comp);
            if(pCheckEdge == NULL)
            {
                // Create the new edges
                ScaffoldEdge* pYZ = new ScaffoldEdge(pZ, dir_yz, comp, dist, 0, 0, SET_INFERRED);
                ScaffoldEdge* pZY = new ScaffoldEdge(pY, dir_zy, comp, dist, 0, 0, SET_INFERRED);
                pYZ->setTwin(pZY);
                pZY->setTwin(pYZ);
                std::cout << "Creating edge: " << *pYZ << " " << *pZY << "\n";
                std::cout << "From: " << *pXY << " " << *pXZ << "\n";
                pY->addEdge(pYZ);
                pZ->addEdge(pZY);
            }
            else
            {
                std::cout << "Edge already exists: " << *pCheckEdge << " inferred dist: " << dist << "\n";
            }

            // Remove the edges pXZ and pZX
            ScaffoldEdge* pZX = pXZ->getTwin();
            pVertex->deleteEdge(pXZ); // deallocates
            pZ->deleteEdge(pZX); // deallocates
            
            pXZ = 0;
            pZX = 0;
            changed_graph = true;
        }
    }

    return changed_graph;
}

void ScaffoldChainVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{

}
