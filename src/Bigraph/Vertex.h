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
#include <set>
#include <vector>
#include <list>
#include <ostream>
#include <iostream>
#include <iterator>
#include "GraphCommon.h"

// Forward declare
class Edge;

// Default edge sorting function, by ID
struct EdgeIDComp
{
	bool operator()(const Edge* pA, const Edge* pB);
};

// Typedefs
typedef std::map<EdgeDesc, Edge*> EdgePtrMap;
typedef std::vector<Edge*> EdgePtrVec;
typedef std::set<EdgeDesc> EdgeDescSet;
typedef std::list<Edge*> EdgePtrList;
typedef EdgePtrMap::iterator EdgePtrMapIter;
typedef EdgePtrMap::const_iterator EdgePtrMapConstIter;
typedef EdgePtrVec::iterator EdgePtrVecIter;
typedef EdgePtrList::iterator EdgePtrListIter;
typedef EdgePtrList::const_iterator EdgePtrListConstIter;

class Vertex
{
	public:
	
		Vertex(VertexID id) : m_id(id), m_color(GC_WHITE) {}
		virtual ~Vertex();

		// Edge list operations
		void addEdge(Edge* ep);
		void removeEdge(Edge* pEdge);
		void removeEdge(const EdgeDesc& ed);
		void deleteEdges();
		void sweepEdges(GraphColor c);
		bool hasEdge(Edge* pEdge) const;
		bool hasEdge(const EdgeDesc& ed) const;
		Edge* getEdge(const EdgeDesc& ed);
		
		// Virtual functions
		virtual void merge(Edge* pEdge);
		virtual void validate() const;
		virtual void sortAdjList();
		
		// getters
		EdgePtrListIter findEdge(const EdgeDesc& ed);
		EdgePtrListConstIter findEdge(const EdgeDesc& ed) const;

		EdgePtrVec findEdgesTo(VertexID id);
		EdgePtrVec getEdges(EdgeDir dir);
		EdgePtrVec getEdges();
		size_t countEdges() const;
		size_t countEdges(EdgeDir dir);

		// 
		VertexID getID() const { return m_id; }
		GraphColor getColor() const { return m_color; }
		void setColor(GraphColor c) { m_color = c; }
		void setEdgeColors(GraphColor c);

		//
		virtual size_t getMemSize() const
		{
			return sizeof(*this) + (m_edges.size() * sizeof(Edge*));
		}

		// Output edges in graphviz format
		void writeEdges(std::ostream& out, int dotFlags) const;

	protected:

		VertexID m_id;
		EdgePtrList m_edges;
		GraphColor m_color;
};

#endif
