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

// Enums

// Arbitrary colors that can be used to indicate the state of a vertex
enum VertexColor
{
	VC_WHITE,
	VC_GRAY,
	VC_BLACK
};

// Typedefs
typedef set<Edge> EdgeSet;
typedef vector<Edge> EdgeVec;
typedef EdgeSet::iterator EdgeSetIter;
typedef EdgeVec::iterator EdgeVecIter;

class Vertex
{
	public:
		Vertex(VertexID id) : m_id(id), m_color(VC_WHITE) {}
		virtual ~Vertex();

		// Add an edge
		void addEdge(Edge e);

		// Add edges in a set
		void addEdges(const EdgeVec& ev);
		
		// Remove an edge
		void removeEdge(Edge e);

		// Check for the precense of an edge
		bool hasEdge(Edge e) const;

		// Merge the data of another vertex into this vertex
		virtual void merge(const Vertex* pV2, const Edge& e);

		// Get the cost of travelling through this node
		virtual int cost() const { return 1; }

		// Set the color of the vertex
		void setColor(VertexColor c) { m_color = c; }

		// Get the color
		VertexColor getColor() const { return m_color; }

		// Find edges to the specified vertex
		EdgeVec findEdgesTo(VertexID id) const;

		// Get the edges in a particular direction
		EdgeVec getEdges(EdgeDir dir) const;

		// Get the edges
		EdgeVec getEdges() const;

		// Count the edges
		size_t countEdges() const { return m_edges.size(); }
		size_t countEdges(EdgeDir dir) const;

		// Return the vert's id
		VertexID getID() const { return m_id; }

		// Output
		friend ostream& operator<<(std::ostream& out, const Vertex& obj);

		// Output edges in graphviz format
		void writeEdges(ostream& out) const;

	private:

		EdgeVec m_mergeRec;
		VertexID m_id;
		EdgeSet m_edges;
		VertexColor m_color;
};

#endif
