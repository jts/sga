#ifndef SEQGRAPH_H
#define SEQGRAPH_H

#include <string>
#include <stdio.h>
#include <vector>
#include "Vertex.h"

using namespace std;

typedef vector<Vertex*> VertexPtrVec;

class SeqGraph
{
	public:
		SeqGraph();
		~SeqGraph();

		// Add a vertex
		void addVertex(Vertex* pVert);
		
		// Remove a vertex
		void removeVertex(VertexID id);

		// Get a vertex
		Vertex* getVertex(VertexID id);

		// Add an edge
		void addEdge(VertexID id1, VertexID id2, EdgeDir dir, EdgeComp comp);

		// Remove an edge
		void removeEdge(VertexID id1, VertexID id2, EdgeDir dir, EdgeComp comp);

		// Dump the graph to a dot file
		void writeDot(string filename);


	private:

		VertexPtrVec m_vertices;
};

#endif
