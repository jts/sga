#include "IVertex.h"

//
//
//
IVertex::~IVertex()
{

}

//
// Add an edge
//
void IVertex::addEdge(VertexID ep, EdgeDir dir, EdgeComp comp)
{
	// Construct the edge
	Edge e(ep, dir, comp);
	m_edges.insert(e);
}

//
// Remove an edge
//
void IVertex::removeEdge(VertexID ep, EdgeDir dir, EdgeComp comp)
{
	// Create the edge that should be removed
	Edge e(ep, dir, comp);
	if(m_edges.find(e) == m_edges.end())
	{
		cerr << "removeEdge:: edge not found " << e << endl;
	}
	m_edges.erase(e);
}

//
// Output
// 
ostream& operator<<(std::ostream& out, const IVertex& obj)
{
	out << obj.m_id << " Edges: \n";
	copy(obj.m_edges.begin(), obj.m_edges.end(), ostream_iterator<Edge>(out, "\n"));
	return out;
}

//
// Output edges in graphviz format
//
void IVertex::writeEdges(ostream& out) const
{
	EdgeSet::const_iterator iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		string color = (iter->getEdgeDir() == ED_SENSE) ? "black" : "red";
		out << "\"" << getID() << "\" -> \"" << iter->getEndpoint();
		out << "\" [color=\"" << color << "\"];\n";
	}
}

