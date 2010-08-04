//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldVertex - A node in a scaffold graph. 
//
#include "ScaffoldVertex.h"

//
ScaffoldVertex::ScaffoldVertex(VertexID id, size_t seqLen) : m_id(id), m_seqLen(seqLen)
{

}

//
ScaffoldVertex::~ScaffoldVertex()
{
    for(ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
          iter != m_edges.end(); ++iter)
    {
        delete *iter;
        *iter = 0;
    }
}

//
void ScaffoldVertex::addEdge(ScaffoldEdge* pEdge)
{
    m_edges.push_back(pEdge);
}

//
ScaffoldEdge* ScaffoldVertex::findEdgeTo(VertexID id, ScaffoldEdgeType type)
{
    for(ScaffoldEdgePtrVector::iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == id && (*iter)->getType() == type)
            return *iter;
    }

    return NULL;
}

//
VertexID ScaffoldVertex::getID() const
{
    return m_id;
}

//
size_t ScaffoldVertex::getNumEdges() const
{
    return m_edges.size();
}

//
size_t ScaffoldVertex::getSeqLen() const
{
    return m_seqLen;
}

//
void ScaffoldVertex::writeEdgesDot(std::ostream* pWriter) const
{
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        *pWriter << "\"" << (*iter)->getStartID() << "\" -> \"" << (*iter)->getEndID();
        std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
        std::string label = ((*iter)->getComp() == EC_SAME) ? "S" : "F";
        *pWriter << "\" [";
        *pWriter << "label=\"" << (*iter)->getDistance() << "\" ";
        *pWriter << "color=\"" << color << "\" ";
        *pWriter << "];";
        *pWriter << "\n";
    }
}
