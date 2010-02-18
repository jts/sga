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

TransitiveGroup::TransitiveGroup(Vertex* pVertex, Edge* pIrreducible) : m_pVertex(pVertex), m_pIrreducibleEdge(pIrreducible)
{

}

//
void TransitiveGroup::addTransitive(Edge* pTransitive)
{
	m_transitivePVec.push_back(pTransitive);
}

//
size_t TransitiveGroup::numTransitive() const
{
	return m_transitivePVec.size();
}

//
Edge* TransitiveGroup::getIrreducible() const
{
	return m_pIrreducibleEdge;
}

//
Edge* TransitiveGroup::getTransitive(size_t idx) const
{
	assert(idx < m_transitivePVec.size());
	return m_transitivePVec[idx];
}

//
Edge* TransitiveGroup::getTransitive(EdgeDesc ed) const
{
	for(size_t i = 0; i < m_transitivePVec.size(); ++i)
	{
		if(m_transitivePVec[i]->getDesc() == ed)
			return m_transitivePVec[i];
	}

	return NULL;
}

// Return the MultiOverlap corresponding to this transitive group
MultiOverlap TransitiveGroup::getMultiOverlap() const
{
	MultiOverlap mo(m_pVertex->getID(), m_pVertex->getSeq());
	mo.add(m_pIrreducibleEdge->getEnd()->getSeq(), m_pIrreducibleEdge->getOverlap());

	for(size_t i = 0; i < m_transitivePVec.size(); ++i)
	{
		Edge* pEdge = m_transitivePVec[i];
		mo.add(pEdge->getEnd()->getSeq(), pEdge->getOverlap());
	}
	return mo;
}
//
bool TransitiveGroup::hasTransitive(EdgeDesc ed) const
{
	Edge* pe = getTransitive(ed);
	return pe != NULL;
}

//
void TransitiveGroup::print() const
{
	std::cout << "\tIrr: " << *m_pIrreducibleEdge << "\n";
	for(size_t i = 0; i < m_transitivePVec.size(); ++i)
	{
		std::cout << "\tTrans: " << *m_transitivePVec[i] << "\n";
	}
}

//
void TransitiveGroup::printMultiOverlap() const
{
	MultiOverlap mo = getMultiOverlap();
	mo.print();
}
