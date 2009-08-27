//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#ifndef VERTEX_H
#define VERTEX_H

// Includes
#include <stdio.h>
#include <map>
#include <vector>
#include <ostream>
#include <iostream>
#include <iterator>
#include "GraphCommon.h"

// Forward declare
class Edge;

// Enums
enum VertexColor
{
	VC_WHITE,
	VC_GRAY,
	VC_BLACK
};

// Typedefs
typedef std::map<EdgeDesc, Edge*> EdgePtrMap;
typedef std::vector<Edge*> EdgePtrVec;
typedef EdgePtrMap::iterator EdgePtrMapIter;
typedef EdgePtrMap::const_iterator EdgePtrMapConstIter;
typedef EdgePtrVec::iterator EdgePtrVecIter;

class Vertex
{
	public:
	
		Vertex(VertexID id) : m_id(id), m_color(VC_WHITE) {}
		virtual ~Vertex();

		// Add an edge
		void addEdge(Edge* ep);
		void removeEdge(Edge* pEdge);
		void removeEdge(const EdgeDesc& ed);
		void deleteEdges();
		bool hasEdge(Edge* pEdge) const;
		bool hasEdge(const EdgeDesc& ed) const;
		Edge* getEdge(const EdgeDesc& ed);
	
		virtual void merge(const Edge* pEdge);
		virtual void validate() const;

		EdgePtrVec findEdgesTo(VertexID id);
		EdgePtrVec getEdges(EdgeDir dir);
		EdgePtrVec getEdges();
		size_t countEdges() const;
		size_t countEdges(EdgeDir dir);

		// 
		VertexID getID() const { return m_id; }
		VertexColor getColor() const { return m_color; }
		void setColor(VertexColor c) { m_color = c; }

		// Output edges in graphviz format
		void writeEdges(std::ostream& out, int dotFlags) const;

	protected:

		VertexID m_id;
		EdgePtrMap m_edges;
		VertexColor m_color;
};

#endif
