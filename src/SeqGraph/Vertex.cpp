#include "Vertex.h"

//
//
//
Vertex::~Vertex()
{

}

//
// Add an edge
//
void Vertex::addEdge(VertexID ep, EdgeDir dir, EdgeComp comp)
{
	// Construct the edge
	Edge e(getID(), ep, dir, comp);
	addEdge(e);
}

//
// Add an edge
//
void Vertex::addEdge(Edge e)
{
	m_edges.insert(e);
}

//
// Add edges from a set
//
void Vertex::addEdges(const EdgeVec& ev)
{
	m_edges.insert(ev.begin(), ev.end());
}

//
// Remove an edge
//
void Vertex::removeEdge(Edge e)
{
	// Check if the edge exists
	if(!hasEdge(e))
	{
		cerr << "removeEdge:: edge not found " << e << endl;
	}
	m_edges.erase(e);
}

//
// Merge data
//
void Vertex::merge(const Vertex* /*pV2*/, const Edge& e)
{
	m_mergeRec.push_back(e);
}

//
// Check for the precense of a particular edge
//
bool Vertex::hasEdge(Edge e) const
{
	return m_edges.find(e) != m_edges.end();
}

//
// Find edges to the specified vertex
//
EdgeVec Vertex::findEdgesTo(VertexID id) const
{
	EdgeSet::const_iterator iter = m_edges.begin();
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
// Find edges in a particular direction
//
EdgeVec Vertex::getEdges(EdgeDir dir) const
{
	EdgeSet::const_iterator iter = m_edges.begin();
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

//
// Get all the edges (as a vector)
//
EdgeVec Vertex::getEdges() const
{
	EdgeVec ev(m_edges.begin(), m_edges.end());
	return ev;
}

//
// Count the edges in a particular direction
// 
size_t Vertex::countEdges(EdgeDir dir) const
{
	EdgeVec ev = getEdges(dir);
	return ev.size();
}

//
// Output
// 
ostream& operator<<(std::ostream& out, const Vertex& obj)
{
	out << obj.m_id << " Edges: \n";
	copy(obj.m_edges.begin(), obj.m_edges.end(), ostream_iterator<Edge>(out, "\n"));
	return out;
}

//
// Output edges in graphviz format
//
void Vertex::writeEdges(ostream& out) const
{
	EdgeSet::const_iterator iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		string color = (iter->getDir() == ED_SENSE) ? "black" : "red";
		string label = (iter->getComp() == EC_NATURAL) ? "0" : "1";
		out << "\"" << iter->getStart() << "\" -> \"" << iter->getEnd();
		out << "\" [color=\"" << color << "\" ";
		out << "label=\"" << label << "\"];\n";
	}
}

