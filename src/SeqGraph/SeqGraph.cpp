#include <assert.h>
#include <ostream>
#include <fstream>
#include <iostream>
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
Vertex* SeqGraph::getVertex(VertexID id) const
{
	VertexPtrMap::const_iterator iter = m_vertices.find(id);
	if(iter == m_vertices.end())
	{
		std::cerr << "Cannot find vertex " << id << " aborting\n";
		assert(iter != m_vertices.end());
	}
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
// Remove an edge
//
void SeqGraph::removeEdge(const Edge& e)
{
	Vertex* pVert1 = getVertex(e.getStart());
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
// Get the IDs of the vertices that do not branch (both sense/antisense degree <= 1)
//
VertexIDVec SeqGraph::getNonBranchingVertices() const
{
	VertexIDVec out;
	VertexPtrMap::const_iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		int senseEdges = iter->second->countEdges(ED_SENSE);
		int antisenseEdges = iter->second->countEdges(ED_ANTISENSE);
		if(antisenseEdges <= 1 && senseEdges <= 1)
		{
			out.push_back(iter->second->getID());
		}
	}
	return out;
}


//
// Get all the paths corresponding to the linear components of the graph
// Precondition: all vertices have in/out degree at most 1 (no branching)
//
PathVector SeqGraph::getLinearComponents()
{
	PathVector outPaths;
	setColors(VC_WHITE);
	VertexPtrMap::iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		// Output the linear path containing this vertex if it hasnt been visited already
		if(iter->second->getColor() != VC_BLACK)
		{
			outPaths.push_back(constructLinearPath(iter->second->getID()));
		}
	}
	assert(checkColors(VC_BLACK));
	return outPaths;
}

//
// Return all the path of nodes that can be linearally reached from this node
// The path expands in both directions so the first node in the path is not necessarily the source
//
Path SeqGraph::constructLinearPath(VertexID id)
{
	Path sensePath;
	Path antisensePath;
	followLinear(id, ED_SENSE, sensePath);
	followLinear(id, ED_ANTISENSE, antisensePath);

	// Construct the final path 
	Path final = reversePath(antisensePath);
	final.insert(final.end(), sensePath.begin(), sensePath.end());
	return final;
}

//
// Recursively follow the graph in the specified direction without branching
// outPath is an out-parameter of the edges that were followed
//
void SeqGraph::followLinear(VertexID id, EdgeDir dir, Path& outPath)
{
	Vertex* pVertex = getVertex(id);
	EdgeVec edges = pVertex->getEdges(dir);

	// Color the vertex
	pVertex->setColor(VC_BLACK);
	
	if(edges.size() == 1)
	{
		Edge& se = edges.front();
		outPath.push_back(se);
		// Correct the direction for the comp
		assert(se.getDir() == dir);
		EdgeDir corrected_dir = correctDir(se.getDir(), se.getComp());

		// Recurse
		followLinear(se.getEnd(), corrected_dir, outPath);
	}
}

//
// Visit each vertex and call the function pointer
// The function returns TRUE if it modifies the graph
//
bool SeqGraph::visit(VertexVisitFunction f)
{
	bool modified = false;
	VertexPtrMap::const_iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		modified = f(this, iter->second) || modified;
	}
	return modified;
}

//
// Set all the vertices in the graph to the given color
//
void SeqGraph::setColors(VertexColor c)
{
	VertexPtrMap::iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		iter->second->setColor(c);
	}
}

//
// Check if all the vertices in the graph are the given color
//
bool SeqGraph::checkColors(VertexColor c)
{
	VertexPtrMap::iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		if(iter->second->getColor() != c)
			return false;
	}
	return true;
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
	std::ofstream out(filename.c_str());
	
	out << "digraph G\n{\n";
	VertexPtrMap::const_iterator iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		VertexID id = iter->second->getID();
		out << id << " [ label =\"" << id << "\"];\n";
		iter->second->writeEdges(out);
	}
	out << "}\n";
	out.close();
}

//
// Reverse the given path
//
Path reversePath(const Path& p)
{
	Path out;
	for(Path::const_reverse_iterator iter = p.rbegin(); iter != p.rend(); ++iter)
		out.push_back(iter->getTwin());

	return out;
}
