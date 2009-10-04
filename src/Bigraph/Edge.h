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
#include "Util.h"
#include "GraphCommon.h"

// Foward declare
class Vertex;

class Edge
{
	public:
		Edge(Vertex* start, Vertex* end, EdgeDir dir, EdgeComp comp) : 
				m_pStart(start), m_pEnd(end), m_dir(dir), m_comp(comp), m_pTwin(NULL), m_color(GC_WHITE) {}

		virtual ~Edge() {}

		void setStart(Vertex* pVert)
		{
			m_pStart = pVert;
		}

		void setTwin(Edge* pEdge)
		{
			m_pTwin = pEdge;
		}

		void setColor(GraphColor c)
		{
			m_color = c;
		}

		Edge* getTwin() const
		{
			assert(m_pTwin != NULL);
			return m_pTwin;
		}
		
		// Merge the given edge with this edge
		virtual void merge(const Edge* pEdge);

		EdgeDesc getTwinDesc() const;

		// Make the direction of the edge that is in the same direction as the current edge
		// but originating in the endpoint vertex
		// This is the transitive direction start --> end *-->* 
		EdgeDir getTransitiveDir() const
		{
			return (m_comp == EC_SAME) ? m_dir : !m_dir;
		}

		// Make the direction of the edge that its twin should point along 
		// start   --->   end
		//       * <--- *
		EdgeDir getTwinDir() const
		{
			return (m_comp == EC_SAME) ? !m_dir : m_dir;
		}

		EdgeDesc getDesc() const
		{
			return EdgeDesc(getEndID(), getDir(), getComp());
		}
		

		// Flip the edge
		void flipComp() { m_comp = !m_comp; }
		void flipDir() { m_dir = !m_dir; }

		virtual void flip()
		{
			flipComp();
			flipDir();
		}

		// Getters
		Vertex* getStart() const;
		Vertex* getEnd() const;
		VertexID getStartID() const; 
		VertexID getEndID() const; 
		EdgeDir getDir() const; 
		EdgeComp getComp() const; 
		bool isSelf() const;
		GraphColor getColor() const;
		virtual std::string getLabel() const { return std::string(""); }

		// Output
		friend std::ostream& operator<<(std::ostream& out, const Edge& obj);

	private:
		
		Edge() {}; // Default constructor is not allowed

		Vertex* m_pStart;
		Vertex* m_pEnd;
		EdgeDir m_dir;
		EdgeComp m_comp;
		Edge* m_pTwin;
		GraphColor m_color;
};

#endif
