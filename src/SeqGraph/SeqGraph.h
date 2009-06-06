#ifndef SEQGRAPH_H
#define SEQGRAPH_H

#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include "Vertex.h"

using namespace std;

// 
// Forward declare
//
class SeqGraph;

//
// Typedefs
//
typedef map<VertexID, Vertex*> VertexPtrMap;
typedef VertexPtrMap::iterator VertexPtrMapIter;
typedef bool(*VertexVisitFunction)(SeqGraph*, Vertex*);
typedef std::vector<Edge> Path;
typedef std::vector<Path> PathVector;
typedef std::vector<VertexID> VertexIDVec;

//
// Functions
//
Path reversePath(const Path& p);

class SeqGraph
{
	public:
		SeqGraph();
		~SeqGraph();

		// Add a vertex
		void addVertex(Vertex* pVert);
		
		// Remove a vertex
		void removeVertex(VertexID id);

		// Check if a vertex exists
		bool hasVertex(VertexID id);

		// Get a vertex
		Vertex* getVertex(VertexID id) const;

		// Add an edge
		void addEdge(const Edge& e);

		// Remove an edge
		void removeEdge(const Edge& e);

		// Merge vertices
		void mergeVertices(VertexID id1, VertexID id2);

		// Simplify the graph by removing transitive edges
		void simplify();

		// Validate that the graph is sane
		void validate();

		// Flip a given vertex
		void flip(VertexID id);

		// Get the IDs of the vertices that do not branch (both sense/antisense degree <= 1)
		VertexIDVec getNonBranchingVertices() const;

		// Get the linear components of a non-branching graph
		PathVector getLinearComponents();

		// Return all the path of nodes that can be linearally reached from this node
		// The path expands in both directions so the first node in the path is not necessarily the source
		Path constructLinearPath(VertexID id);

		// Print simple summary statistics to stdout
		void stats() const;
		
		// Visit each vertex in the graph and perform the visit function
		bool visit(VertexVisitFunction f);

		// Dump the graph to a dot file
		void writeDot(string filename) const;

	private:

		// Set/check the colors for the entire graph
		void setColors(VertexColor c);
		bool checkColors(VertexColor c);

		void followLinear(VertexID id, EdgeDir dir, Path& outPath);

		// Merge two vertices along a specified edge
		void mergeAlongEdge(Vertex* pV1, Vertex* pV2, const Edge& edge);

		// Vertex collection
		VertexPtrMap m_vertices;
};

#endif
