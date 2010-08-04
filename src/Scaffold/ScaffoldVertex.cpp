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
VertexID ScaffoldVertex::getID() const
{
    return m_id;
}

