#ifndef EDGE_H
#define EDGE_H

#include <ostream>
#include "Util.h"
#include "Common.h"

using namespace std;

class Edge
{
	public:
		Edge() {};
		Edge(VertexID start, VertexID end, EdgeDir dir, EdgeComp comp, int overlap) : 
				m_start(start), m_end(end), m_dir(dir), m_comp(comp), m_overlap(overlap) {}

		// Generate the twin edge, the edge from the end node back to the start node
		Edge getTwin() const
		{
			return Edge(m_end, m_start, getTwinDir(), m_comp, m_overlap);
		}

		// Make the direction of the edge that its twin should point along 
		EdgeDir getTwinDir() const
		{
			return (m_comp == EC_SAME) ? !m_dir : m_dir;
		}

		// Flip the edge
		void flip()
		{
			m_comp = !m_comp;
			m_dir = !m_dir;
		}

		// Getters
		VertexID getStart() const { return m_start; }
		VertexID getEnd() const { return m_end; }
		EdgeDir getDir() const { return m_dir; }
		EdgeComp getComp() const { return m_comp; }
		int getOverlap() const { return m_overlap; }
		bool isSelf() const { return m_start == m_end; }
		double getWeight() const { return m_weight; }

		// Setters
		void setWeight(double w) { m_weight = w; }

		// Equality operator
		bool operator==(const Edge& obj) const
		{
			return (m_start == obj.m_start) && (m_end == obj.m_end) && 
					(m_dir == obj.m_dir) && (m_comp == obj.m_comp);
		}

		// Less than
		bool operator<(const Edge& obj) const
		{
			if(m_start < obj.m_start)
				return true;
			else if(m_start > obj.m_start)
				return false;
			else if(m_end < obj.m_end)
				return true;
			else if(m_end > obj.m_end)
				return false;
			else if(m_dir < obj.m_dir)
				return true;
			else if(m_dir > obj.m_dir)
				return false;
			else if(m_comp < obj.m_comp)
				return true;
			else if(m_comp > obj.m_comp)
				return false;
			return false;
		}

		// Output
		friend ostream& operator<<(std::ostream& out, const Edge& obj)
		{
			out << obj.m_start << "," << obj.m_end << "," << obj.m_dir << "," << obj.m_comp;
			return out;
		}

	private:
		VertexID m_start;
		VertexID m_end;
		EdgeDir m_dir;
		EdgeComp m_comp;
		int m_overlap;
		double m_weight;

};

#endif
