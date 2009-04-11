#ifndef IVERTEX_H
#define IVERTEX_H

// Includes
#include <stdio.h>
#include <set>
#include <vector>
#include <ostream>
#include <iostream>
#include <iterator>
#include "Common.h"
#include "Edge.h"

// Namespace
using namespace std;

// Typedefs
typedef set<Edge> EdgeSet;
typedef vector<Edge> EdgeVec;
typedef EdgeSet::iterator EdgeSetIter;
typedef EdgeVec::iterator EdgeVecIter;

class Vertex
{
	public:
		Vertex(VertexID id) : m_id(id) {}
		virtual ~Vertex();

		// Add an edge
		void addEdge(VertexID ep, EdgeDir dir, EdgeComp comp);
		void addEdge(Edge e);

		// Add edges in a set
		void addEdges(const EdgeVec& ev);
		
		// Remove an edge
		void removeEdge(Edge e);

		// Check for the precense of an edge
		bool hasEdge(Edge e) const;

		// Find edges to the specified vertex
		EdgeVec findEdgesTo(VertexID id) const;

		// Get the edges in a particular direction
		EdgeVec getEdgesInDir(EdgeDir dir) const;

		// Get the edges
		EdgeVec getEdges() const;

		// Count the edges in a particular direction
		size_t countEdgesInDir(EdgeDir dir) const;

		// Return the vert's id
		VertexID getID() const { return m_id; }

		// Output
		friend ostream& operator<<(std::ostream& out, const Vertex& obj);

		// Output edges in graphviz format
		void writeEdges(ostream& out) const;

	private:

		VertexID m_id;
		EdgeSet m_edges;
};

#endif
