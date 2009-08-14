#ifndef VERTEX_H
#define VERTEX_H

// Includes
#include <stdio.h>
#include <set>
#include <vector>
#include <ostream>
#include <iostream>
#include <iterator>
#include "GraphCommon.h"
#include "Edge.h"

// Namespace
using namespace std;

// Enums

// Arbitrary colors that can be used to indicate the state of a vertex
enum VertexColor
{
	VC_WHITE,
	VC_GRAY,
	VC_BLACK
};


template<typename VD, typename ET>
class Vertex
{
	public:

		// Typedefs
		typedef VD VertexData;
		typedef ET EdgeType;
		typedef set<EdgeType> EdgeSet;
		typedef vector<EdgeType> EdgeVec;
		typedef typename EdgeSet::iterator EdgeSetIter;
		typedef typename EdgeSet::iterator EdgeSetConstIter;
		typedef typename EdgeVec::iterator EdgeVecIter;
	
		Vertex(VertexID id) : m_id(id), m_color(VC_WHITE) {}
		Vertex(VertexID id, VD d) : m_id(id), m_color(VC_WHITE), m_data(d) {}
		~Vertex() {};

		//
		// Add an edge
		//
		void addEdge(EdgeType e)
		{
			std::pair<EdgeSetIter, bool> result = m_edges.insert(e);
			if(!result.second)
			{
				//std::cerr << "Warning, added duplicate edge " << e << std::endl;
			}
		}

		//
		// Add edges in a set
		//
		void addEdges(const EdgeVec& ev)
		{
			m_edges.insert(ev.begin(), ev.end());
		}
		
		//
		// Remove an edge
		//
		void removeEdge(EdgeType e)
		{
			// Check if the edge exists
			if(!hasEdge(e))
			{
				cerr << "removeEdge:: edge not found " << e << " in vertex " << *this <<  endl;
				assert(false);
			}
			m_edges.erase(e);
		}

		//
		// Check for the precense of an edge
		//
		bool hasEdge(EdgeType e) const
		{
			return m_edges.find(e) != m_edges.end();
		}

		//
		// Return the matching edge
		//
		EdgeType getEdge(EdgeType e) const
		{
			 EdgeSetConstIter i = m_edges.find(e);
			 assert(i != m_edges.end());
			 return *i;
		}

		//
		// Get the cost of travelling through this node
		//
		int cost() const { return 1; }

		//
		// Set the color of the vertex
		//
		void setColor(VertexColor c) { m_color = c; }

		//
		// Get the color
		//
		VertexColor getColor() const { return m_color; }

		//
		// Find edges to the specified vertex
		//
		EdgeVec findEdgesTo(VertexID id) const
		{
			EdgeSetConstIter iter = m_edges.begin();
			EdgeVec outEdges;
			for(; iter != m_edges.end(); ++iter)
			{
				if(iter->getEnd() == id)
				{
					outEdges.push_back(*iter);
				}
			}
			return outEdges;
		}


		//
		// Get the edges in a particular direction
		//
		EdgeVec getEdges(EdgeDir dir) const
		{
			EdgeSetConstIter iter = m_edges.begin();
			EdgeVec outEdges;
			for(; iter != m_edges.end(); ++iter)
			{
				if(iter->getDir() == dir)
				{
					outEdges.push_back(*iter);
				}
			}
			return outEdges;
		}
		

		// Get the edges
		EdgeVec getEdges() const
		{
			EdgeVec ev(m_edges.begin(), m_edges.end());
			return ev;
		}


		// Count the edges
		size_t countEdges() const { return m_edges.size(); }
		size_t countEdges(EdgeDir dir) const
		{
			EdgeVec ev = getEdges(dir);
			return ev.size();
		}
		

		// Return the vert's id
		VertexID getID() const { return m_id; }

		// Return the data this vertex carries
		VD& getData() { return m_data; }

		// Output
		friend ostream& operator<<(std::ostream& out, const Vertex& obj)
		{
			out << obj.m_id << " Edges: \n";
			copy(obj.m_edges.begin(), obj.m_edges.end(), ostream_iterator<EdgeType>(out, "\n"));
			return out;
		}	


		// Output edges in graphviz format
		void writeEdges(ostream& out, int dotFlags) const
		{
			EdgeSetConstIter iter = m_edges.begin();
			for(; iter != m_edges.end(); ++iter)
			{
				if(dotFlags & DF_UNDIRECTED)
				{
					if(iter->getStart() < iter->getEnd())
					{
						out << "\"" << iter->getStart() << "\" -- \"" << iter->getEnd() << "\"";
					}
				}
				else
				{
					out << "\"" << iter->getStart() << "\" -> \"" << iter->getEnd();
					string color = (iter->getDir() == ED_SENSE) ? "black" : "red";
					string label = (iter->getComp() == EC_SAME) ? "S" : "F";
					out << "\" [color=\"" << color << "\" ";
					out << "label=\"" << label << " (" << iter->getData() << ") \"];";
				}
				out << "\n";

			}
		}
		

	private:

		EdgeVec m_mergeRec;
		VertexID m_id;
		EdgeSet m_edges;
		VertexColor m_color;
		VD m_data;
};

#endif
