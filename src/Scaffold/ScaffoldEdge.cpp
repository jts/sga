//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldEdge - An edge in a scaffold graph. 
//
#include "ScaffoldEdge.h"
#include "ScaffoldVertex.h"

ScaffoldEdge::ScaffoldEdge(ScaffoldVertex* pEnd, 
                           EdgeDir dir, EdgeComp comp, 
                           int distance, double stdDev, 
                           int numPairs, ScaffoldEdgeType type) : m_pEnd(pEnd), m_pTwin(NULL),
                                                                  m_distance(distance), m_stdDev(stdDev),
                                                                  m_numPairs(numPairs), m_type(type)
{
    m_edgeData.setDir(dir);
    m_edgeData.setComp(comp);
}

//
void ScaffoldEdge::setTwin(ScaffoldEdge* pTwin)
{
    m_pTwin = pTwin;
}

//
VertexID ScaffoldEdge::getStartID() const
{
    assert(m_pTwin != NULL);
    return m_pTwin->getEndID();
}

//
VertexID ScaffoldEdge::getEndID() const
{
    return m_pEnd->getID();
}

//
EdgeDir ScaffoldEdge::getDir() const
{
    return m_edgeData.getDir();
}

//
EdgeComp ScaffoldEdge::getComp() const
{
    return m_edgeData.getComp();
}

//
int ScaffoldEdge::getDistance() const
{
    return m_distance;
}

//
double ScaffoldEdge::getStdDev() const
{
    return m_stdDev;
}

//
ScaffoldEdgeType ScaffoldEdge::getType() const
{
    return m_type;
}
