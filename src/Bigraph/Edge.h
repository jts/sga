//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Base bidirectional edge class 
//

#ifndef EDGE_H
#define EDGE_H

#include <ostream>
#include <boost/pool/object_pool.hpp>
#include "Match.h"
#include "Util.h"
#include "GraphCommon.h"
#include "Vertex.h"

class Edge
{
	public:
		Edge(Vertex* start, Vertex* end, EdgeDir dir, EdgeComp comp, SeqCoord m) : 
				m_pStart(start), m_pEnd(end), 
				m_pTwin(NULL), m_matchCoord(m), 
				m_dir(dir), m_comp(comp), m_color(GC_WHITE) {}

		~Edge() {}
		
		// High level modification functions
		
		// Join merges pEdge into this edge, with pEdge describing the starting point
		void join(const Edge* pEdge);

		// Extend merged pEdge into this edge, with pEdge describing the endpoint
		void extend(const Edge* pEdge);

		// Post merge update function
		void update() {}

		// Sequence Coordinate functions
		
		// update functions
		void extendMatchFullLength();
		void extendMatch(int ext_len);
		void offsetMatch(int offset);
		void updateSeqLen(int newLen);

		// access functions
		size_t getSeqLen() const;
		size_t getMatchLength() const { return m_matchCoord.length(); }
		const SeqCoord& getMatchCoord() const { return m_matchCoord; }
		std::string getMatchStr() const;
		Match getMatch() const;		

		// setters
		void setStart(Vertex* pVert) { 	m_pStart = pVert; }
		void setTwin(Edge* pEdge) { m_pTwin = pEdge; }
		void setColor(GraphColor c) { m_color = c; }

		// getters
		VertexID getStartID() const { return m_pStart->getID(); }
		VertexID getEndID() const { return m_pEnd->getID(); }
		inline Vertex* getStart() const { return m_pStart; }
		inline Vertex* getEnd() const { return m_pEnd; }
		inline EdgeDir getDir() const { return m_dir; }
 		inline EdgeComp getComp() const { return m_comp; }		
		inline Edge* getTwin() const { assert(m_pTwin != NULL); return m_pTwin; }
		EdgeDesc getTwinDesc() const;
		std::string getLabel() const;

		inline bool isSelf() const { return m_pStart == m_pEnd; }
		inline GraphColor getColor() const { return m_color; }
		size_t getMemSize() const { return sizeof(*this); }

		// Returns the direction of an edge that continues in the same direction
		// as this edge, corrected for complementary 
		inline EdgeDir getTransitiveDir() const { return (m_comp == EC_SAME) ? m_dir : !m_dir; }

		// Make the direction of the edge that its twin should point along 
		inline EdgeDir getTwinDir() const { return (m_comp == EC_SAME) ? !m_dir : m_dir; }
		inline EdgeDesc getDesc() const { return EdgeDesc(getEndID(), getDir(), getComp()); }
		
		// Flip the edge
		inline void flipComp() { m_comp = !m_comp; }
		inline void flipDir() { m_dir = !m_dir; }
		void flip() { flipComp(); flipDir(); }

		// Memory management
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

		// Validate that the edge is sane
		void validate() const;

		// Output
		friend std::ostream& operator<<(std::ostream& out, const Edge& obj);

	protected:
			
		// Return the pointer to the memory pool
		// The first time this function is called the memory pool is created
		static boost::object_pool<Edge>* mempool() 
		{
			if(!m_spMempool)
			{
				m_spMempool = new boost::object_pool<Edge>;
			}
			return m_spMempool;
    	}

  		static boost::object_pool<Edge>* m_spMempool; 

		Edge() {}; // Default constructor is not allowed

		Vertex* m_pStart;
		Vertex* m_pEnd;
		Edge* m_pTwin;
		SeqCoord m_matchCoord;
		EdgeDir m_dir;
		EdgeComp m_comp;
		GraphColor m_color;
};

#endif
