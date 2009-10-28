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
	EdgePtrListIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		delete *iter;
		*iter = NULL;
	}
}

// Add an edge
void Vertex::addEdge(Edge* ep)
{
	assert(ep->getStart() == this);

#ifdef VALIDATE
	for(EdgePtrListConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		if((*iter)->getEndID() == ep->getEndID())
		{
			std::cout << "Attempted to add duplicate edge with ID: " << ep->getEndID() << " to vertex: " << ep->getStartID() << "\n";
			std::cout << "Added in desc: " << ep->getDesc() << " curr desc: " << (*iter)->getDesc() << "\n";
			assert(false);
		}
	}
#endif
	m_edges.push_back(ep);
}

// Remove an edge
void Vertex::removeEdge(Edge* pEdge)
{
	// Check if the edge exists
	assert(pEdge != NULL);
	removeEdge(pEdge->getDesc());
}

// Remove edge
// Note - this does not delete the edge through the pointer
void Vertex::removeEdge(const EdgeDesc& ed)
{

	EdgePtrListIter iter = findEdge(ed);
	assert(iter != m_edges.end());
	m_edges.erase(iter);
}

// Delete all the edges, and their twins, from this vertex
void Vertex::deleteEdges()
{
	EdgePtrListIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		Edge* pEdge = *iter;
		Edge* pTwin = pEdge->getTwin();
		Vertex* pPartner = pEdge->getEnd();
		pPartner->removeEdge(pTwin);
		delete pEdge;
		pEdge = NULL;
		delete pTwin;
		pTwin = NULL;
		*iter = NULL;
	}
	m_edges.clear();
}

// Delete edges that are marked
// This only deletes the edge and not its twin
void Vertex::sweepEdges(GraphColor c)
{
	EdgePtrListIter iter = m_edges.begin();
	while(iter != m_edges.end())
	{
		Edge* pEdge = *iter;
		if(pEdge->getColor() == c)
		{
			delete pEdge;
			pEdge = NULL;
			iter = m_edges.erase(iter);
		}
		else
			++iter;
	}
}

// Return the iterator to the edge matching edgedesc
EdgePtrListIter Vertex::findEdge(const EdgeDesc& ed)
{
	for(EdgePtrListIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		if((*iter)->getDesc() == ed)
			return iter;
	}
	return m_edges.end();
}

//
EdgePtrListConstIter Vertex::findEdge(const EdgeDesc& ed) const
{
	for(EdgePtrListConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		if((*iter)->getDesc() == ed)
			return iter;
	}
	return m_edges.end();
}

//
void Vertex::sortAdjList()
{
	EdgeIDComp comp;
	m_edges.sort(comp);
}

//
void Vertex::validate() const
{
	EdgeDescSet edSet;
	// Ensure the twin edge exists for every edge
	for(EdgePtrListConstIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
	{
		Edge* pEdge = *iter;
		Edge* pTwinEdge = pEdge->getTwin();
		Vertex* pEndpoint = pEdge->getEnd();
		if(pTwinEdge == NULL)
			std::cerr << "Warning, twin pointer for edge " << *pEdge << " is NULL\n";
		else if(!pEndpoint->hasEdge(pEdge->getTwinDesc()))
			std::cerr << "Warning edge " << *pEdge << " does not have a twin with desc " << pEdge->getTwinDesc() << "\n";
		
		std::pair<EdgeDescSet::iterator, bool> result = edSet.insert(pEdge->getDesc());
		if(result.second != true)
		{
			std::cerr << "Error: edge with description " << pEdge->getDesc() 
			          << " is in the adj list twice\n";
			assert(false);
		}
	}
}

// Merge
void Vertex::merge(Edge* /*pEdge*/)
{
	// base case does nothing
}

// Check for the presence of an edge
bool Vertex::hasEdge(Edge* pEdge) const
{
	return hasEdge(pEdge->getDesc());
}

//
bool Vertex::hasEdge(const EdgeDesc& ed) const
{
	return findEdge(ed) != m_edges.end();
}

// Return the edge matching the descriptions
Edge* Vertex::getEdge(const EdgeDesc& ed)
{
	 EdgePtrListIter i = findEdge(ed);
	 assert(i != m_edges.end());
	 return *i;
}

// Find edges to the specified vertex
EdgePtrVec Vertex::findEdgesTo(VertexID id)
{
	EdgePtrListConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	for(; iter != m_edges.end(); ++iter)
	{
		if((*iter)->getEndID() == id)
			outEdges.push_back(*iter);
	}
	return outEdges;
}


//
// Get the edges in a particular direction
// This preserves the ordering of the edges
//
EdgePtrVec Vertex::getEdges(EdgeDir dir)
{
	EdgePtrListConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	for(; iter != m_edges.end(); ++iter)
	{
		if((*iter)->getDir() == dir)
			outEdges.push_back(*iter);
	}
	return outEdges;
}


// Get the edges
EdgePtrVec Vertex::getEdges()
{
	EdgePtrListConstIter iter = m_edges.begin();
	EdgePtrVec outEdges;
	for(; iter != m_edges.end(); ++iter)
		outEdges.push_back(*iter);
	return outEdges;	
}

void Vertex::setEdgeColors(GraphColor c) 
{ 
	for(EdgePtrListIter iter = m_edges.begin(); iter != m_edges.end(); ++iter)
		(*iter)->setColor(c);
}

// Count the edges
// This function is not necessarily constant time
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
	EdgePtrListConstIter iter = m_edges.begin();
	for(; iter != m_edges.end(); ++iter)
	{
		if(dotFlags & DF_UNDIRECTED)
		{
			if((*iter)->getStartID() < (*iter)->getEndID())
			{
				out << "\"" << (*iter)->getStart() << "\" -- \"" << (*iter)->getEnd() << "\"";
			}
		}
		else
		{
			out << "\"" << (*iter)->getStartID() << "\" -> \"" << (*iter)->getEndID();
			std::string color = ((*iter)->getDir() == ED_SENSE) ? "black" : "red";
			std::string label = ((*iter)->getComp() == EC_SAME) ? "S" : "F";
			out << "\" [color=\"" << color << "\" ";
			out << "label=\"" << (*iter)->getLabel() << "\"];";
		}
		out << "\n";
	}
}

bool EdgeIDComp::operator()(const Edge* pA, const Edge* pB) 
{
   	return pA->getEndID() < pB->getEndID();
}
