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

// Edge sorting function, by length
struct EdgeLenComp
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
typedef EdgePtrVec::const_iterator EdgePtrVecConstIter;
typedef EdgePtrList::iterator EdgePtrListIter;
typedef EdgePtrList::const_iterator EdgePtrListConstIter;

class Vertex
{
	public:
	
		Vertex(VertexID id, const std::string& s) : m_id(id), 
												   	m_seq(s), 
													m_readCount(1), 
													m_pPairVertex(NULL),
													m_color(GC_WHITE) {}
		~Vertex();

		// High-level modification functions
		
		// Merge another vertex into this vertex, as specified by pEdge
		void merge(Edge* pEdge);

		// sort the edges by the ID of the vertex they point to
		void sortAdjListByID();

		// sort the edges by the length of the label of the edge
		void sortAdjListByLen();

		// Ensure that all the edges are unique
		void makeUnique(); 

		// Edge list operations
		void addEdge(Edge* ep);
		void removeEdge(Edge* pEdge);
		void removeEdge(const EdgeDesc& ed);
		void deleteEdges();
		void sweepEdges(GraphColor c);
		bool hasEdge(Edge* pEdge) const;
		bool hasEdge(const EdgeDesc& ed) const;

		Edge* getEdge(const EdgeDesc& ed);
		EdgePtrVec findEdgesTo(VertexID id);
		EdgePtrVec getEdges(EdgeDir dir);
		EdgePtrVec getEdges();
		EdgePtrVecIter findEdge(const EdgeDesc& ed);
		EdgePtrVecConstIter findEdge(const EdgeDesc& ed) const;

		size_t countEdges() const;
		size_t countEdges(EdgeDir dir);

		// Ensure the vertex data is sane
		void validate() const;
		
		// setters
		void setPairVertex(Vertex* pPair);
		void clearPairVertex();
		void setColor(GraphColor c) { m_color = c; }
		void setEdgeColors(GraphColor c);

		// getters
		VertexID getID() const { return m_id; }
		GraphColor getColor() const { return m_color; }
		Vertex* getPairVertex() const { return m_pPairVertex; }
		size_t getReadCount() const { return m_readCount; }
		const std::string& getSeq() const { return m_seq; }
		size_t getMemSize() const;

		// Output edges in graphviz format
		void writeEdges(std::ostream& out, int dotFlags) const;

	private:

		// Ensure all the edges in DIR are unique
		void makeUnique(EdgeDir dir, EdgePtrVec& uniqueVec);

		VertexID m_id;
		EdgePtrVec m_edges;
		std::string m_seq;
		int m_readCount;
		Vertex* m_pPairVertex;
		GraphColor m_color;
};

#endif
