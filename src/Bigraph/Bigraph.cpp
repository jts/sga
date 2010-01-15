//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Bidirectional graph 
//
#include <assert.h>
#include <ostream>
#include <fstream>
#include <iostream>
#include "Bigraph.h"
#include "Timer.h"
#include <malloc.h>

//
//
//
Bigraph::Bigraph()
{
}

//
//
//
Bigraph::~Bigraph()
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
void Bigraph::addVertex(Vertex* pVert)
{
	m_vertices.insert(std::make_pair(pVert->getID(), pVert));
}

//
// Remove a vertex that is not connected to any other
//
void Bigraph::removeIslandVertex(Vertex* pVertex)
{
	assert(pVertex->countEdges() == 0);

	// Remove the vertex from the collection
	VertexID id = pVertex->getID();
	delete pVertex;
	m_vertices.erase(id);
}

//
// Remove a vertex
//
void Bigraph::removeConnectedVertex(Vertex* pVertex)
{
	// Remove the edges pointing to this Vertex
	pVertex->deleteEdges();

	// Remove the vertex from the collection
	VertexID id = pVertex->getID();
	delete pVertex;
	m_vertices.erase(id);
}


//
// Check for the existance of a vertex
//
bool Bigraph::hasVertex(VertexID id)
{
	VertexPtrMapIter iter = m_vertices.find(id);
	return iter != m_vertices.end();
}

//
// Get a vertex
//
Vertex* Bigraph::getVertex(VertexID id) const
{
	VertexPtrMapConstIter iter = m_vertices.find(id);
	if(iter == m_vertices.end())
		return NULL;
	return iter->second;
}

//
// Add an edge
//
void Bigraph::addEdge(Edge* pEdge)
{
	Vertex* pVert = pEdge->getStart();
	pVert->addEdge(pEdge);
}

//
// Remove an edge
//
void Bigraph::removeEdge(const EdgeDesc& ed)
{
	Vertex* pVert = getVertex(ed.id);
	pVert->removeEdge(ed);
}

//
// High level merge function that does not specify an edge
//
void Bigraph::mergeVertices(VertexID id1, VertexID id2)
{
	Vertex* pVert1 = getVertex(id1);

	// Get the edges from vertex1 to vertex2
	EdgePtrVec edgesTo = pVert1->findEdgesTo(id2);

	if(edgesTo.empty())
	{
		std::cerr << "mergeVertices: vertices are not connected\n";
		return;
	}

	if(edgesTo.size() > 1)
	{
		std::cerr << "mergeVertces: cannot merge because of ambigious edges\n";
		return;
	}

	// There is a single unique edge between the vertices
	Edge* mergeEdge = *edgesTo.begin();

	// Call the real merging function
	merge(mergeEdge);
}

//
// Merge two vertices along the specified edge
//
void Bigraph::merge(Edge* pEdge)
{
	Vertex* pV1 = pEdge->getStart();
	Vertex* pV2 = pEdge->getEnd();
	//std::cout << "Merging " << pV1->getID() << " with " << pV2->getID() << "\n";

	// Merge the data
	pV1->merge(pEdge);

	// Get the twin edge (the edge in v2 that points to v1)
	Edge* pTwin = pEdge->getTwin();

	// Ensure v2 has the twin edge
	assert(pV2->hasEdge(pTwin));

	// Get the edge set opposite of the twin edge (which will be the new edges in this direction for V1)
	EdgePtrVec transEdges = pV2->getEdges(!pTwin->getDir());

	// Move the edges from pV2 to pV1
	for(EdgePtrVecIter iter = transEdges.begin(); iter != transEdges.end(); ++iter)
	{
		Edge* pTransEdge = *iter;

		// Remove the edge from V2, this does not destroy the edge
		pV2->removeEdge(pTransEdge);

		// Join pEdge to the start of transEdge
		// This updates the starting point of pTransEdge to be V1
		// This calls Edge::extend on the twin edge
		pTransEdge->join(pEdge);
		assert(pTransEdge->getDir() == pEdge->getDir());
		pV1->addEdge(pTransEdge); // add to V1

		// Notify the edges they have been updated
		pTransEdge->update();
		pTransEdge->getTwin()->update();
	}

	// Remove the edge from pV1 to pV2
	pV1->removeEdge(pEdge);
	delete pEdge;
	pEdge = 0;

	// Remove the edge from pV2 to pV1
	pV2->removeEdge(pTwin);
	delete pTwin;
	pEdge = 0;

	// Remove V2
	// It is guarenteed to not be connected
	removeIslandVertex(pV2);
	//validate();
}

//
void Bigraph::sweepVertices(GraphColor c)
{
	VertexPtrMapIter iter = m_vertices.begin();
	while(iter != m_vertices.end())
	{
		VertexPtrMapIter next = iter;
		++next;
		if(iter->second->getColor() == c)
			removeConnectedVertex(iter->second);
		iter = next;
	}
}


//
void Bigraph::sweepEdges(GraphColor c)
{
	for(VertexPtrMapIter iter = m_vertices.begin(); iter != m_vertices.end(); ++iter)
		iter->second->sweepEdges(c);
}

//	Simplify the graph by compacting singular edges
void Bigraph::simplify()
{
	simplify(ED_SENSE);
	simplify(ED_ANTISENSE);
}

// Simplify the graph by compacting edges in the given direction
void Bigraph::simplify(EdgeDir dir)
{
	Timer* pTimer = new Timer("simplify timer");
	bool graph_changed = true;
	while(graph_changed)
	{
		int proc_count = 0;
		graph_changed = false;
		VertexPtrMapIter iter = m_vertices.begin(); 
		while(iter != m_vertices.end())
		{
			// Get the edges for this direction
			EdgePtrVec edges = iter->second->getEdges(dir);

			// If there is a single edge in this direction, merge the vertices
			// Don't merge singular self edges though
			if(edges.size() == 1 && !edges.front()->isSelf())
			{
				// Check that the edge back is singular as well
				Edge* pSingle = edges.front();
				Edge* pTwin = pSingle->getTwin();
				Vertex* pV2 = pSingle->getEnd();
				if(pV2->countEdges(pTwin->getDir()) == 1)
				{
					merge(pSingle);
					graph_changed = true;
				}
			}

			if(proc_count++ % 20000 == 0)
			{
				printf("processed_count\t%d\tsimplify_time\t%lf\n", proc_count, pTimer->getElapsedTime());
				pTimer->reset();
			}
			++iter;
		}
	} 
	delete pTimer;
}

//
// Sort the adjacency list for each vertex
//
void Bigraph::sortVertexAdjLists()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
		iter->second->sortAdjList();
}

//
// Validate that the graph is sane
//
void Bigraph::validate()
{
	VertexPtrMapIter iter = m_vertices.begin();
	for(; iter != m_vertices.end(); ++iter)
	{
		iter->second->validate();
	}
}

//
// Flip a vertex
//
void Bigraph::flip(VertexID /*id*/)
{
	assert(false);
#if 0
	// TODO: update this code
	Vertex* pVertex = getVertex(id);
	EdgePtrVec edges = pVertex->getEdges();

	for(EdgePtrVecIter iter = edges.begin(); iter != edges.end(); ++iter)
	{
		// Get the old twin
		GraphEdgeType twin = iter->getTwin();
		
		GraphEdgeType flipped = *iter; 
		flipped.flip();

		// Remove the edge from the source ver
		pVertex->removeEdge(*iter);
		pVertex->addEdge(flipped);

		// Update the partner by deleting the old twin and 
		Vertex* pV2 = getVertex(twin.getStart());
		pV2->removeEdge(twin);
		pV2->addEdge(flipped.getTwin());
	}
#endif
}

//
// Get the IDs of the vertices that do not branch (both sense/antisense degree <= 1)
//
VertexIDVec Bigraph::getNonBranchingVertices() const
{
	VertexIDVec out;
	VertexPtrMapConstIter iter = m_vertices.begin(); 
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
PathVector Bigraph::getLinearComponents()
{
	PathVector outPaths;
	setColors(GC_WHITE);
	VertexPtrMapIter iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		// Output the linear path containing this vertex if it hasnt been visited already
		if(iter->second->getColor() != GC_BLACK)
		{
			outPaths.push_back(constructLinearPath(iter->second->getID()));
		}
	}
	assert(checkColors(GC_BLACK));
	return outPaths;
}

//
// Return all the path of nodes that can be linearally reached from this node
// The path expands in both directions so the first node in the path is not necessarily the source
//
Path Bigraph::constructLinearPath(VertexID id)
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
void Bigraph::followLinear(VertexID id, EdgeDir dir, Path& outPath)
{
	Vertex* pVertex = getVertex(id);
	EdgePtrVec edges = pVertex->getEdges(dir);

	// Color the vertex
	pVertex->setColor(GC_BLACK);
	
	if(edges.size() == 1)
	{
		Edge* pSingle = edges.front();
		outPath.push_back(pSingle);
		// Correct the direction for the comp
		assert(pSingle->getDir() == dir);
		EdgeDir corrected_dir = correctDir(pSingle->getDir(), pSingle->getComp());

		// Recurse
		followLinear(pSingle->getEndID(), corrected_dir, outPath);
	}
}

//
Path Bigraph::reversePath(const Path& path)
{
	Path out;
    for(Path::const_reverse_iterator iter = path.rbegin(); iter != path.rend(); ++iter)
		out.push_back((*iter)->getTwin());
	return out;
}

//
// Visit each vertex and call the function pointer
// The function returns TRUE if it modifies the graph
//
bool Bigraph::visit(VertexVisitFunction f)
{
	bool modified = false;
	VertexPtrMapConstIter iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		modified = f(this, iter->second) || modified;
	}
	return modified;
}

//
// Set all the vertices and edges in the graph to the given color
//
void Bigraph::setColors(GraphColor c)
{
	VertexPtrMapIter iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		iter->second->setColor(c);
		iter->second->setEdgeColors(c);
	}
}

//
// Check if all the vertices in the graph are the given color
//
bool Bigraph::checkColors(GraphColor c)
{
	VertexPtrMapIter iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		if(iter->second->getColor() != c)
		{
			std::cerr << "Warning vertex " << iter->second->getID() << " is color " << iter->second->getColor() << " expected " << c << "\n";
			return false;
		}
	}
	return true;
}

//
// Output simple stats
//
void Bigraph::stats() const
{
	int numVerts = 0;
	int numEdges = 0;

	VertexPtrMapConstIter iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		numEdges += iter->second->countEdges();
		++numVerts;
	}

	std::cout << "Graph has " << numVerts << " vertices and " << numEdges << " edges\n";
}

//
// Output mem stats
//
void Bigraph::printMemSize() const
{
	size_t numVerts = 0;
	size_t vertMem = 0;

	size_t numEdges = 0;
	size_t edgeMem = 0;

	VertexPtrMapConstIter iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		++numVerts;
		vertMem += iter->second->getMemSize();

		EdgePtrVec edges = iter->second->getEdges();
		for(EdgePtrVecIter edgeIter = edges.begin(); edgeIter != edges.end(); ++edgeIter)
		{
			++numEdges;
			edgeMem += (*edgeIter)->getMemSize();
		}
	}
	printf("num verts: %zu using %zu bytes (%.2lf per vert)\n", numVerts, vertMem, double(vertMem) / numVerts);
	printf("num edges: %zu using %zu bytes (%.2lf per edge)\n", numEdges, edgeMem, double(edgeMem) / numEdges);
	printf("total: %zu\n", edgeMem + vertMem);
}

//
// Dump the graph to a dot file
//
void Bigraph::writeDot(std::string filename, int dotFlags) const
{
	std::ofstream out(filename.c_str());
	
	std::string graphType = (dotFlags & DF_UNDIRECTED) ? "graph" : "digraph";

	out << graphType << " G\n{\n";
	VertexPtrMapConstIter iter = m_vertices.begin(); 
	for(; iter != m_vertices.end(); ++iter)
	{
		VertexID id = iter->second->getID();
		out << "\"" << id << "\" [ label =\"" << id << "\" ";
		out << "];\n";
		iter->second->writeEdges(out, dotFlags);
	}
	out << "}\n";
	out.close();
}

