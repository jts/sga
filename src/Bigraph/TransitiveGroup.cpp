//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// TransitiveGroup - A set of edges that form a 
// mutually transitive set
//
#include "Edge.h"
#include "TransitiveGroup.h"

TransitiveGroup::TransitiveGroup(Vertex* pVertex, Edge* pIrreducible) : m_pVertex(pVertex)
{
	m_edges.push_back(pIrreducible);
}

//
void TransitiveGroup::add(Edge* pEdge)
{
	m_edges.push_back(pEdge);
}

//
size_t TransitiveGroup::numElements() const
{
	return m_edges.size();
}

//
Edge* TransitiveGroup::getIrreducible() const
{
	assert(!m_edges.empty());
	// The irreducible edge is the first one in the group
	return m_edges[0];
}

//
Edge* TransitiveGroup::getEdge(size_t idx) const
{
	assert(idx < m_edges.size());
	return m_edges[idx];
}

//
Edge* TransitiveGroup::getEdge(EdgeDesc ed) const
{
	for(size_t i = 0; i < m_edges.size(); ++i)
	{
		if(m_edges[i]->getDesc() == ed)
			return m_edges[i];
	}

	return NULL;
}

// Return the MultiOverlap corresponding to this transitive group
MultiOverlap TransitiveGroup::getMultiOverlap() const
{
	MultiOverlap mo(m_pVertex->getID(), m_pVertex->getSeq());

	for(size_t i = 0; i < m_edges.size(); ++i)
	{
		Edge* pEdge = m_edges[i];
		mo.add(pEdge->getEnd()->getSeq(), pEdge->getOverlap());
	}
	return mo;
}
//
bool TransitiveGroup::hasEdge(EdgeDesc ed) const
{
	Edge* pe = getEdge(ed);
	return pe != NULL;
}

//
void TransitiveGroup::print() const
{
	std::cout << "\tIrr: " << *m_edges[0] << "\n";
	for(size_t i = 1; i < m_edges.size(); ++i)
	{
		std::cout << "\tTrans: " << *m_edges[i] << "\n";
	}
}

//
void TransitiveGroup::printMultiOverlap() const
{
	MultiOverlap mo = getMultiOverlap();
	mo.printPileup();
}
