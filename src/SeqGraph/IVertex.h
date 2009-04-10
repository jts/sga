#ifndef IVERTEX_H
#define IVERTEX_H

// Includes
#include <stdio.h>
#include <set>
#include <ostream>
#include <iostream>
#include <iterator>
#include "Common.h"
#include "Edge.h"

// Namespace
using namespace std;

// Typedefs
typedef set<Edge> EdgeSet;

class IVertex
{
	public:
		IVertex(VertexID id) : m_id(id) {}
		virtual ~IVertex();

		// Add an edge
		void addEdge(VertexID ep, EdgeDir dir, EdgeComp comp);
		
		// Remove an edge
		void removeEdge(VertexID ep, EdgeDir dir, EdgeComp comp);

		// Return the vert's id
		VertexID getID() const { return m_id; }

		// Output
		friend ostream& operator<<(std::ostream& out, const IVertex& obj);

		// Output edges in graphviz format
		void writeEdges(ostream& out) const;

	private:

		VertexID m_id;
		EdgeSet m_edges;
};

#endif
