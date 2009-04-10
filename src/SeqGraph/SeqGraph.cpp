#include <assert.h>
#include "SeqGraph.h"

//
//
//
SeqGraph::SeqGraph()
{
}

//
//
//
SeqGraph::~SeqGraph()
{
}

//
// Add a vertex
//
void SeqGraph::addVertex(IVertex* pVert)
{
	// Make sure the id is correct
	assert(pVert->getID() == m_vertices.size());
	m_vertices.push_back(pVert);
}

//
// Remove a vertex
//
void SeqGraph::removeVertex(VertexID id)
{
	assert(id < m_vertices.size());
	m_vertices[id] = NULL;
}

//
// Get a vertex
//
IVertex* SeqGraph::getVertex(VertexID id)
{
	assert(id < m_vertices.size());
	IVertex* pVert = m_vertices[id];
	assert(pVert != NULL);
	return pVert;
}

//
// Add an edge
//
void SeqGraph::addEdge(VertexID id1, VertexID id2, EdgeDir dir, EdgeComp comp)
{
	IVertex* pVert1 = getVertex(id1);
	pVert1->addEdge(id2, dir, comp);
}

//
// Remove an edge
//
void SeqGraph::removeEdge(VertexID id1, VertexID id2, EdgeDir dir, EdgeComp comp)
{
	IVertex* pVert1 = getVertex(id1);
	pVert1->removeEdge(id2, dir, comp);

}

//
// Dump the graph to a dot file
//
void SeqGraph::writeDot(string filename)
{
	(void)filename;
	cout << "digraph G\n{\n";
	VertexPtrVec::const_iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		VertexID id = (*iter)->getID();
		cout << id << " [ label =\"" << id << "\"];\n";
		(*iter)->writeEdges(cout);
	}
	cout << "}\n";
}

