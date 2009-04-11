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
void SeqGraph::addVertex(Vertex* pVert)
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

	// Remove the edges pointing to this Vertex
	Vertex* pVertex = getVertex(id);
	EdgeVec ev = pVertex->getEdges();
	for(EdgeVecIter iter = ev.begin(); iter != ev.end(); ++iter)
	{
		Edge twin = iter->getTwin();
		Vertex* pV2 = getVertex(twin.getStart());
		pV2->removeEdge(twin);
	}

	// Remove the vertex from the collection
	delete pVertex;
	m_vertices[id] = NULL;
}

//
// Get a vertex
//
Vertex* SeqGraph::getVertex(VertexID id)
{
	assert(id < m_vertices.size());
	Vertex* pVert = m_vertices[id];
	assert(pVert != NULL);
	return pVert;
}

//
// Add an edge
//
void SeqGraph::addEdge(VertexID id1, VertexID id2, EdgeDir dir, EdgeComp comp)
{
	Vertex* pVert1 = getVertex(id1);
	pVert1->addEdge(id2, dir, comp);
}

//
// Remove an edge
//
void SeqGraph::removeEdge(VertexID id1, VertexID id2, EdgeDir dir, EdgeComp comp)
{
	Vertex* pVert1 = getVertex(id1);
	Edge e(id1, id2, dir, comp);
	pVert1->removeEdge(e);
}

//
// High level merge function that does not specify an edge
//
void SeqGraph::mergeVertices(VertexID id1, VertexID id2)
{
	Vertex* pVert1 = getVertex(id1);
	Vertex* pVert2 = getVertex(id2); 

	// Get the edges from vertex1 to vertex2
	EdgeVec edgesTo = pVert1->findEdgesTo(id2);

	if(edgesTo.empty())
	{
		cerr << "mergeVertices: vertices are not connected\n";
		return;
	}

	if(edgesTo.size() > 1)
	{
		cerr << "mergeVertces: cannot merge because of ambigious edges\n";
		return;
	}

	// There is a single unique edge between the vertices
	Edge mergeEdge = *edgesTo.begin();

	// Call the real merging function
	mergeAlongEdge(pVert1, pVert2, mergeEdge);
}

//
// Merge two vertices along the specified edge
//
void SeqGraph::mergeAlongEdge(Vertex* pV1, Vertex* pV2, const Edge& edge)
{
	// Construct the twin edge (the edge in v2 that points to v1)
	Edge twinEdge = edge.getTwin();

	// Ensure v2 has the twin edge
	assert(pV2->hasEdge(twinEdge));

	// Get the edge set opposite of the twin edge
	EdgeVec transEdges = pV2->getEdgesInDir(!twinEdge.getDir());

	// Should the edges be flipped?
	bool doFlip = (edge.getComp() == EC_REVERSE);

	// Add the new edges to V1
	for(EdgeVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
	{
		// If the verts dont have the same comp, flip the edge
		if(doFlip)
			iter->flip();

		assert(iter->getDir() == edge.getDir());

		// Build the new edge and add it to V1
		Edge e(pV1->getID(), iter->getEnd(), iter->getDir(), iter->getComp());
		pV1->addEdge(e);

		// Add the twin edge to the new partner node
		Edge twin = e.getTwin();
		Vertex* pV3 = getVertex(twin.getStart());
		pV3->addEdge(twin);
	}

	// Remove the edge from pV1 to pV2
	pV1->removeEdge(edge);

	// Remove the edge from pV2 to pV1
	pV2->removeEdge(twinEdge);

	// Check if V2 should be completely deleted
	size_t edgeCount = pV2->countEdgesInDir(twinEdge.getDir());
	if(edgeCount == 0)
	{
		removeVertex(pV2->getID());
	}
}

//
// Dump the graph to a dot file
//
void SeqGraph::writeDot(string filename) const
{
	(void)filename;
	cout << "digraph G\n{\n";
	VertexPtrVec::const_iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		// If the edge has been deleted, skip
		if(*iter == NULL)
			continue;

		VertexID id = (*iter)->getID();
		cout << id << " [ label =\"" << id << "\"];\n";
		(*iter)->writeEdges(cout);
	}
	cout << "}\n";
}

