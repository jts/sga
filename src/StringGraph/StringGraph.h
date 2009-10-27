//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// String Graph - Bidirectional graph of sequence reads
// and their overlaps
// Inherits from Bigraph::Vertex/Bigraph::Edge
//
#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H

#include "Match.h"
#include "Bigraph.h"
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <string>

#define SE_CAST(x) static_cast<StringEdge*>((x))
#define CSE_CAST(x)  static_cast<const StringEdge*>((x))
#define SV_CAST(x) static_cast<StringVertex*>((x))
#define CSV_CAST(x)  static_cast<const StringVertex*>((x))


typedef Bigraph StringGraph;
class StringVertex;
class StringEdge;

// String edge sorting function, by length
struct StringEdgeLenComp
{
	bool operator()(const Edge* pA, const Edge* pB);
};


// Derived from Bigraph Edge
class StringEdge : public Edge
{
	public:

		// constructors
		StringEdge(Vertex* start, 
                   Vertex* end, 
		    	   EdgeDir dir, 
			       EdgeComp comp, 
			       SeqCoord m,
				   int nd) : Edge(start, end, dir, comp), m_matchCoord(m), m_numDiff(nd) {}
		
		// functions
		virtual void flip();
		virtual void join(const Edge* pEdge);
		virtual void extend(const Edge* pEdge);

		// Get a match structure that describes the mapping from V1 to V2
		Match getMatch() const;		
		
		// Match coordinate bookkeeping
		void extendMatch(int ext_len);
		void offsetMatch(int offset);
		void completeMatch();
		void updateSeqLen(int newLen);

		// getters
		virtual std::string getLabel() const;
		size_t getSeqLen() const;
		
		size_t getMatchLength() const { return m_matchCoord.length(); }
		std::string getMatchStr() const;
		const SeqCoord& getMatchCoord() const { return m_matchCoord; }

		int getNumDiff() const { return m_numDiff; }

		// Validate that the edge is sane
		void validate() const;
		
	private:
		
		// The coords of the starting vertex, which match the ending vertex
		SeqCoord m_matchCoord;

		// The number of differences in the matching region
		int m_numDiff;
};

// Derived from Bigraph vertex
class StringVertex : public Vertex
{
	public:
		// constructors
		StringVertex(VertexID id, const std::string& s) : Vertex(id), m_seq(s), m_readCount(1) {}
		
		// functions
		virtual void merge(Edge* pEdge);
		virtual void validate() const;
		virtual void sortAdjList();

		size_t getReadCount() const { return m_readCount; }
		const std::string& getSeq() const { return m_seq; }

	private:
		std::string m_seq;
		size_t m_readCount;
};

#endif
