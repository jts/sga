#ifndef EDGE_H
#define EDGE_H

#include <ostream>
#include "Common.h"

using namespace std;

enum EdgeDir
{	
	ED_SENSE,
	ED_ANTISENSE
};

enum EdgeComp
{
	EC_NATURAL,
	EC_REVERSE
};

class Edge
{
	public:
		Edge(VertexID ep, EdgeDir dir, EdgeComp comp) : m_endpoint(ep), m_dir(dir), m_comp(comp) {}


		// Getters
		VertexID getEndpoint() const { return m_endpoint; }
		EdgeDir getEdgeDir() const { return m_dir; }
		EdgeComp getEdgeComp() const { return m_comp; }

		// Equality operator
		bool operator==(const Edge& obj) const
		{
			return (m_endpoint == obj.m_endpoint) && 
					(m_dir == obj.m_dir) && (m_comp == obj.m_comp);
		}

		// Less than
		bool operator<(const Edge& obj) const
		{
			return (m_endpoint < obj.m_endpoint) ||
					(m_dir < obj.m_dir) || (m_comp < obj.m_comp);

		}

		// Output
		friend ostream& operator<<(std::ostream& out, const Edge& obj)
		{
			out << obj.m_endpoint << "," << obj.m_dir << "," << obj.m_comp;
			return out;
		}

	private:
		VertexID m_endpoint;
		EdgeDir m_dir;
		EdgeComp m_comp;


};

#endif
