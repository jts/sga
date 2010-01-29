//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// StringVertex - Derived from Bigraph/Vertex
//
#ifndef STRINGVERTEX_H
#define STRINGVERTEX_H

#include "StringGraph.h"
class StringEdge;

class StringVertex : public Vertex
{
	public:
		// constructors
		StringVertex(VertexID id, const std::string& s) : Vertex(id), m_seq(s), m_readCount(1), m_pPairVertex(NULL) {}
		virtual ~StringVertex(); 

		// functions
		virtual void merge(Edge* pEdge);
		virtual void validate() const;
		virtual void sortAdjList();

		// Ensure that all the edges are unique
		void makeUnique(); 

		void setPairVertex(StringVertex* pPair) { m_pPairVertex = pPair; }
		void clearPairVertex() { m_pPairVertex = NULL; }

		StringVertex* getPairVertex() const { return m_pPairVertex; }
		size_t getReadCount() const { return m_readCount; }
		const std::string& getSeq() const { return m_seq; }

	private:

		// 
		void makeUnique(EdgeDir dir, EdgePtrVec& uniqueVec);

		// data
		std::string m_seq;
		size_t m_readCount;

		// Pointer to the read pair of this sequence, if there is one
		StringVertex* m_pPairVertex;
};

#endif
