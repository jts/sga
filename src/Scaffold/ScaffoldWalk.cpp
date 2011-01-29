//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldWalk - a walk through a scaffold graph
//
#include "ScaffoldWalk.h"
#include "ScaffoldVertex.h"

ScaffoldWalk::ScaffoldWalk(ScaffoldVertex* pStartVertex) : m_pStartVertex(pStartVertex)
{
}

void ScaffoldWalk::addEdge(ScaffoldEdge* pEdge)
{
    m_edges.push_back(pEdge);
}

//
int64_t ScaffoldWalk::getGapSum() const
{
    int64_t sum = 0;
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        sum += (*iter)->getDistance();

    return sum;
}

//
int64_t ScaffoldWalk::getContigLengthSum() const
{
    int64_t sum = m_pStartVertex->getSeqLen();

    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        sum += (*iter)->getEnd()->getSeqLen();
    return sum;
}

//
ScaffoldVertexPtrVector ScaffoldWalk::getVertices() const
{
    ScaffoldVertexPtrVector outVertices;
    outVertices.push_back(m_pStartVertex);
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
        outVertices.push_back((*iter)->getEnd());
    return outVertices;
}


//
int ScaffoldWalk::findVertex(ScaffoldVertex* pVertex) const
{
    if(pVertex == m_pStartVertex)
        return 0;

    int idx = 1;
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
                                              iter != m_edges.end();
                                              ++iter)
    {
        if(pVertex == (*iter)->getEnd())
            return idx;
        ++idx;
    }
    return -1;
}

// Print the edges in dot format
void ScaffoldWalk::printDot(std::ostream& out) const
{
    ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
    for(; iter != m_edges.end(); ++iter)
    {
        ScaffoldEdge* pEdge = *iter;
        out << "\"" << pEdge->getStartID() << "\" -> \"" << pEdge->getEndID() << "\" [color=\"blue\"];\n";
    }
}

void ScaffoldWalk::print() const
{
    for(ScaffoldEdgePtrVector::const_iterator iter = m_edges.begin();
                                              iter != m_edges.end();
                                              ++iter)
    {
        std::cout << *(*iter) << " ";
    }
    std::cout << "\n";
}
