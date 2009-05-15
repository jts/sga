#ifndef SEQGRAPH_H
#define SEQGRAPH_H

#include <string>
#include <stdio.h>
#include <vector>
#include "Vertex.h"

using namespace std;

typedef vector<Vertex*> VertexPtrVec;
typedef VertexPtrVec::iterator VertexPtrVecIter;

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
		void addEdge(const Edge& e);

		// Merge vertices
		void mergeVertices(VertexID id1, VertexID id2);

		// Simplify the graph by removing transitive edges
		void simplify();

		// Validate that the graph is sane
		void validate();

		// Flip a given vertex
		void flip(VertexID id);

		// Print simple summary statistics to stdout
		void stats() const;
		
		// Dump the graph to a dot file
		void writeDot(string filename) const;


	private:

		// Merge two vertices along a specified edge
		void mergeAlongEdge(Vertex* pV1, Vertex* pV2, const Edge& edge);

		// Vertex collection
		VertexPtrVec m_vertices;
};

#endif
