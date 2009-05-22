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
	VertexPtrMap::iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		delete iter->second;
		iter->second = NULL;
	}
}

//
// Add a vertex
//
void SeqGraph::addVertex(Vertex* pVert)
{
	// Make sure the id is correct
	m_vertices.insert(std::make_pair(pVert->getID(), pVert));
}

//
// Remove a vertex
//
void SeqGraph::removeVertex(VertexID id)
{
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
	m_vertices.erase(id);
}

//
// Check for the existance of a vertex
//
bool SeqGraph::hasVertex(VertexID id)
{
	VertexPtrMap::iterator iter = m_vertices.find(id);
	return iter != m_vertices.end();
}

//
// Get a vertex
//
Vertex* SeqGraph::getVertex(VertexID id)
{
	VertexPtrMap::iterator iter = m_vertices.find(id);
	assert(iter != m_vertices.end());
	return iter->second;
}

//
// Add an edge
//
void SeqGraph::addEdge(const Edge& e)
{
	Vertex* pVert1 = getVertex(e.getStart());
	pVert1->addEdge(e);
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
	// Merge the data
	pV1->merge(pV2, edge);

	// Construct the twin edge (the edge in v2 that points to v1)
	Edge twinEdge = edge.getTwin();

	// Ensure v2 has the twin edge
	assert(pV2->hasEdge(twinEdge));

	// Get the edge set opposite of the twin edge
	EdgeVec transEdges = pV2->getEdges(!twinEdge.getDir());

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
		Edge e(pV1->getID(), iter->getEnd(), iter->getDir(), iter->getComp(), iter->getOverlap());
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
	size_t edgeCount = pV2->countEdges(twinEdge.getDir());
	if(edgeCount == 0)
	{
		removeVertex(pV2->getID());
	}
}

//
//	Simplify the graph by removing transitive edges
//
void SeqGraph::simplify()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		for(int d = 0; d < ED_COUNT; ++d)
		{
			EdgeDir dir = EDGE_DIRECTIONS[d];
			
			// Get the edges for this direction
			EdgeVec edges = iter->second->getEdges(dir);

			// If there is a single edge in this direction, merge the vertices
			// Don't merge singular self edges though
			if(edges.size() == 1 && !edges.front().isSelf())
			{
				Edge single = edges.front();
				mergeAlongEdge(iter->second, getVertex(single.getEnd()), single);
			}
		}
	}
}

//
// Validate that the graph is sane
//
void SeqGraph::validate()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		// Ensure the twin edge exists for every edge
		EdgeVec edges = iter->second->getEdges();
		for(EdgeVecIter iter = edges.begin(); iter != edges.end(); ++iter)
		{
			Vertex* pV2 = getVertex(iter->getEnd());
			if(!pV2->hasEdge(iter->getTwin()))
			{
				cerr << "Warning, twin edge does not exist for " << iter->getStart() 
						<< "," << iter->getEnd() << endl;
			}
		}
	}

}

//
// Flip a vertex
//
void SeqGraph::flip(VertexID id)
{
	Vertex* pVertex = getVertex(id);
	EdgeVec edges = pVertex->getEdges();

	for(EdgeVecIter iter = edges.begin(); iter != edges.end(); ++iter)
	{
		// Get the old twin
		Edge twin = iter->getTwin();
		
		Edge flipped = *iter; 
		flipped.flip();

		// Remove the edge from the source ver
		pVertex->removeEdge(*iter);
		pVertex->addEdge(flipped);

		// Update the partner by deleting the old twin and 
		Vertex* pV2 = getVertex(twin.getStart());
		pV2->removeEdge(twin);
		pV2->addEdge(flipped.getTwin());
	}
}

//
// Output simple stats
//
void SeqGraph::stats() const
{
	int numVerts = 0;
	int numEdges = 0;

	VertexPtrMap::const_iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		numEdges += iter->second->countEdges();
		++numVerts;
	}

	std::cout << "Graph has " << numVerts << " vertices and " << numEdges << " edges\n";
}

//
// Dump the graph to a dot file
//
void SeqGraph::writeDot(string filename) const
{
	(void)filename;
	cout << "digraph G\n{\n";
	VertexPtrMap::const_iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		VertexID id = iter->second->getID();
		cout << id << " [ label =\"" << id << "\"];\n";
		iter->second->writeEdges(cout);
	}
	cout << "}\n";
}

