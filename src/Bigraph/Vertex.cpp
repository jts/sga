//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Vertex - Generic vertex class for bigraph
//
#include "Vertex.h"
#include "Edge.h"

Vertex::~Vertex()
{
	EdgePtrMapIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		delete iter->second;
		iter->second = NULL;
	}
}

// Add an edge
void Vertex::addEdge(Edge* ep)
{
	std::pair<EdgePtrMapIter, bool> result = m_edges.insert(std::make_pair(ep->getDesc(), ep));
	if(!result.second)
	{
		std::cerr << "Error attempt to add duplicate edge " << ep << std::endl;
		assert(false);
	}
}

// Remove an edge
void Vertex::removeEdge(Edge* pEdge)
{
	// Check if the edge exists
	assert(pEdge != NULL);
	removeEdge(pEdge->getDesc());
}

// Remove edge
// Note - this does not delete the edge
void Vertex::removeEdge(const EdgeDesc& ed)
{
	EdgePtrMapIter iter = m_edges.find(ed);
	assert(iter != m_edges.end());
	m_edges.erase(iter);
}

// Merge
void Vertex::merge(const Edge* pEdge)
{
	// Signal all the vertices connected to this one that it has been updated
	for(EdgePtrMapIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		iter->second->getEnd()->partnerUpdate(iter->second->getTwin(), pEdge);
	}
}

// Partner update handler
void Vertex::partnerUpdate(const Edge* /*pPartner*/, const Edge* /*pMerged*/)
{
	// do nothing in the base case
}

// Check for the presence of an edge
bool Vertex::hasEdge(Edge* pEdge) const
{
	return hasEdge(pEdge->getDesc());
}

//
bool Vertex::hasEdge(const EdgeDesc& ed) const
{
	return m_edges.find(ed) != m_edges.end();
}

// Return the edge matching the descriptions
Edge* Vertex::getEdge(const EdgeDesc& ed)
{
	 EdgePtrMapConstIter i = m_edges.find(ed);
	 assert(i != m_edges.end());
	 return i->second;
}

// Find edges to the specified vertex
EdgePtrVec Vertex::findEdgesTo(VertexID id)
{
	EdgePtrMapConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	for(; iter != m_edges.end(); ++iter)
	{
		if(iter->second->getEndID() == id)
		{
			outEdges.push_back(iter->second);
		}
	}
	return outEdges;
}


//
// Get the edges in a particular direction
//
EdgePtrVec Vertex::getEdges(EdgeDir dir)
{
	EdgePtrMapConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	outEdges.reserve(m_edges.size());
	for(; iter != m_edges.end(); ++iter)
	{
		if(iter->second->getDir() == dir)
		{
			outEdges.push_back(iter->second);
		}
	}
	return outEdges;
}


// Get the edges
EdgePtrVec Vertex::getEdges()
{
	EdgePtrMapConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	outEdges.reserve(m_edges.size());
	for(; iter != m_edges.end(); ++iter)
		outEdges.push_back(iter->second);
	return outEdges;	
}


// Count the edges
size_t Vertex::countEdges() const
{ 
	return m_edges.size(); 
}

//
size_t Vertex::countEdges(EdgeDir dir)
{
	EdgePtrVec ev = getEdges(dir);
	return ev.size();
}

// Output edges in graphviz format
void Vertex::writeEdges(std::ostream& out, int dotFlags) const
{
	EdgePtrMapConstIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		if(dotFlags & DF_UNDIRECTED)
		{
			if(iter->second->getStartID() < iter->second->getEndID())
			{
				out << "\"" << iter->second->getStart() << "\" -- \"" << iter->second->getEnd() << "\"";
			}
		}
		else
		{
			out << "\"" << iter->second->getStartID() << "\" -> \"" << iter->second->getEndID();
			std::string color = (iter->second->getDir() == ED_SENSE) ? "black" : "red";
			std::string label = (iter->second->getComp() == EC_SAME) ? "S" : "F";
			out << "\" [color=\"" << color << "\" ";
			out << "label=\"" << iter->second->getLabel() << "\"];";
		}
		out << "\n";
	}
}

