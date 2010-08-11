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

ScaffoldEdge::ScaffoldEdge(ScaffoldVertex* pEnd, ScaffoldLink link) : m_pEnd(pEnd), m_pTwin(NULL), m_link(link)
{
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
ScaffoldVertex* ScaffoldEdge::getEnd() const
{
    return m_pEnd;
}

//
ScaffoldEdge* ScaffoldEdge::getTwin() const
{
    return m_pTwin;
}

//
EdgeDir ScaffoldEdge::getDir() const
{
    return m_link.edgeData.getDir();
}

//
EdgeComp ScaffoldEdge::getComp() const
{
    return m_link.edgeData.getComp();
}

//
int ScaffoldEdge::getDistance() const
{
    return m_link.distance;
}

//
double ScaffoldEdge::getStdDev() const
{
    return m_link.stdDev;
}

//
ScaffoldLinkType ScaffoldEdge::getType() const
{
    return m_link.type;
}

char ScaffoldEdge::getTypeCode() const
{
    switch(m_link.type)
    {
        case SLT_DISTANCEEST:
            return 'D';
        case SLT_REFERENCE:
            return 'R';
        case SLT_INFERRED:
            return 'I';
        default:
            return 'N';
    }
}


std::string ScaffoldEdge::makeLinkString() const
{
    std::stringstream ss;
    ss << getEndID() << "," << getDistance() << "," << getStdDev() << "," <<
          getDir() << "," << getComp() << "," << getTypeCode();

    return ss.str();
}

//
std::ostream& operator<<(std::ostream& out, ScaffoldEdge& edge)
{
    out << edge.getStartID() << " -- " << edge.getEndID() << "," << edge.getDistance() << "," << edge.getStdDev() 
        << "," << edge.getDir() << "," << edge.getComp() << "," << edge.getTypeCode();
    return out;
}

//
bool ScaffoldEdgePtrDistanceCompare(ScaffoldEdge* pXY, ScaffoldEdge* pXZ)
{
    return pXY->getDistance() < pXZ->getDistance();
}

