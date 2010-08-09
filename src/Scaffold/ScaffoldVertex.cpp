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
ScaffoldEdge* ScaffoldVertex::findEdgeTo(VertexID id, ScaffoldEdgeType type) const
{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == id && (*iter)->getType() == type)
            return *iter;
    }

    return NULL;
}

//
ScaffoldEdge* ScaffoldVertex::findEdgeTo(VertexID id, EdgeDir dir, EdgeComp comp) const
{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getEndID() == id && (*iter)->getDir() == dir && (*iter)->getComp() == comp)
            return *iter;
    }

    return NULL;
}


//
ScaffoldEdgePtrVector ScaffoldVertex::getEdges(EdgeDir dir)
{
    ScaffoldEdgePtrVector out;
    for(ScaffoldEdgePtrVector::iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        if((*iter)->getDir() == dir)
            out.push_back(*iter);
    }
    return out;
}

void ScaffoldVertex::deleteEdges()
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        delete *iter;
        *iter = NULL;
        ++iter;
    }
    m_edges.clear();
}

void ScaffoldVertex::deleteEdgesAndTwins()
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        ScaffoldEdge* pEdge = *iter;
        pEdge->getEnd()->deleteEdge(pEdge->getTwin());
        delete pEdge;
        *iter = NULL;
        ++iter;
    }
    m_edges.clear();
}

// Remove an edge from the edge vector and delete it
void ScaffoldVertex::deleteEdge(ScaffoldEdge* pEdge)
{
    ScaffoldEdgePtrVector::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        if(*iter == pEdge)
            break;
        ++iter;
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter);
    delete pEdge;
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
bool ScaffoldVertex::isRepeat() const
{
    return m_classification == SVC_REPEAT;
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
        *pWriter << "\"" << (*iter)->getStartID() << "\" -> \"" << (*iter)->getEndID();
        std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
        *pWriter << "\" [";
        *pWriter << "label=\"" << (*iter)->getDistance() << "\" ";
        *pWriter << "color=\"" << color << "\" ";
        *pWriter << "];";
        *pWriter << "\n";
    }
}
