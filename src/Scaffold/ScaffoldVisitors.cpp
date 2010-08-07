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

void ScaffoldAStatisticVisitor::postvisit(ScaffoldGraph* /*pGraph*/)
{
    std::cerr << "A-statistic: classified " << m_numUnique << " components as unique\n";
    std::cerr << "A-statistic: classified " << m_numRepeat << " components as repeat\n"; 
}
