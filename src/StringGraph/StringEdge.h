//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// StringEdge class, derived from Bigraph Edge
//
#ifndef STRINGEDGE_H
#define STRINGEDGE_H

#include "StringGraph.h"
class StringVertex;

// String edge sorting function, by length
struct StringEdgeLenComp
{
	bool operator()(const Edge* pA, const Edge* pB);
};

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
		virtual void update();

		// Get a match structure that describes the mapping from V1 to V2
		Match getMatch() const;		
		
		// Match coordinate bookkeeping
		void extendMatch(int ext_len);
		void offsetMatch(int offset);
		void updateSeqLen(int newLen);

		// Update the number of differences between the sequences. Must be called after a 
		// merge/extend
		void updateDifferenceCount();
		
		// Make the match full length
		void extendMatchFullLength();

		// getters
		virtual std::string getLabel() const;
		size_t getSeqLen() const;
		
		size_t getMatchLength() const { return m_matchCoord.length(); }
		std::string getMatchStr() const;
		const SeqCoord& getMatchCoord() const { return m_matchCoord; }

		int getNumDiff() const { return m_numDiff; }

		// Validate that the edge is sane
		void validate() const;

		// Memory allocation/deallocation, delegate to singleton pool
		void* operator new(size_t /*size*/)
		{
			return mempool()->malloc();
		}
		/*
		void* operator new(size_t size, StringEdge* pEdge)
		{
			return mempool()->malloc();
		}
		*/
  		void operator delete(void* target, size_t /*size*/)
		{
			(void)target;
			//return mempool()->free((StringEdge*)target);
		}

	private:

		// Return the pointer to the memory pool
		// The first time this function is called the memory pool is created
		static boost::object_pool<StringEdge>* mempool() 
		{
			if(!m_spMempool)
			{
				m_spMempool = new boost::object_pool<StringEdge>;
			}
			return m_spMempool;
    	}

  		static boost::object_pool<StringEdge>* m_spMempool; 

		// The coords of the starting vertex, which match the ending vertex
		SeqCoord m_matchCoord;

		// The number of differences in the matching region
		int m_numDiff;
};

#endif
