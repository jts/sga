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
ScaffoldVertex::ScaffoldVertex(VertexID id, size_t seqLen) : m_id(id), m_seqLen(seqLen), m_classification(SVC_UNKNOWN)
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
void ScaffoldVertex::setAStatistic(double v)
{
    m_AStatistic = v;
}

//
void ScaffoldVertex::setClassification(ScaffoldVertexClassification classification)
{
    m_classification = classification;
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
double ScaffoldVertex::getAStatistic() const
{
    return m_AStatistic;
}

//
ScaffoldVertexClassification ScaffoldVertex::getClassification() const
{
    return m_classification;
}

//
std::string ScaffoldVertex::getColorString() const
{
    switch(m_classification)
    {
        case SVC_UNKNOWN:
            return "gray";
        case SVC_UNIQUE:
            return "white";
        case SVC_REPEAT:
            return "red";
        default:
            return "white";
    }
}

//
void ScaffoldVertex::writeDot(std::ostream* pWriter) const
{
   VertexID id = getID();
   *pWriter << "\"" << id << "\" [ label =\"" << id << "," << getSeqLen() << "\" ";
   *pWriter << "style=\"filled\" fillcolor=\"" << getColorString() << "\" ";
   *pWriter << "];\n";
   writeEdgesDot(pWriter);
}

//
void ScaffoldVertex::writeEdgesDot(std::ostream* pWriter) const
{
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        if((*iter)->getStartID() < (*iter)->getEndID())
        {
            *pWriter << "\"" << (*iter)->getStartID() << "\" -- \"" << (*iter)->getEndID();
            std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
            *pWriter << "\" [";
            *pWriter << "label=\"" << (*iter)->getDistance() << "\" ";
            *pWriter << "color=\"" << color << "\" ";
            *pWriter << "];";
            *pWriter << "\n";
        }
    }
}
