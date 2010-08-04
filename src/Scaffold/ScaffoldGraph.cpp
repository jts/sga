//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldGraph - A graph representing long-distance
// relationships between contigs.
//
#include "ScaffoldGraph.h"
#include "SeqReader.h"

ScaffoldGraph::ScaffoldGraph()
{

}

//
ScaffoldGraph::~ScaffoldGraph()
{
    for(ScaffoldVertexMap::iterator iter = m_vertices.begin();
         iter != m_vertices.end(); ++iter)
    {
        delete iter->second;
        iter->second = NULL;
    }
}

//
void ScaffoldGraph::addVertex(ScaffoldVertex* pVertex)
{
    m_vertices.insert(std::make_pair(pVertex->getID(), pVertex));
}

//
void ScaffoldGraph::addEdge(ScaffoldVertex* pVertex, ScaffoldEdge* pEdge)
{
    assert(pVertex != NULL);
    pVertex->addEdge(pEdge);
}

//
ScaffoldVertex* ScaffoldGraph::getVertex(VertexID id) const
{
    ScaffoldVertexMap::const_iterator iter = m_vertices.find(id);
    if(iter == m_vertices.end())
        return NULL;
    return iter->second;
}

// 
void ScaffoldGraph::loadVertices(const std::string& filename)
{
    SeqReader reader(filename);
    SeqRecord sr;
    while(reader.get(sr))
    {
        ScaffoldVertex* pVertex = new ScaffoldVertex(sr.id, sr.seq.length());
        addVertex(pVertex);
    }    
}

//
void ScaffoldGraph::writeDot(const std::string& outFile) const
{
    (void)outFile;
}
