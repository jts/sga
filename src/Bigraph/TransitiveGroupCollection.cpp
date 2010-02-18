//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// TransitiveGroupCollection - A collection
// of TransitiveGroups, representing the partitioning
// of edges for a given vertex
//
#include "TransitiveGroupCollection.h"
#include "Edge.h"

//
TransitiveGroupCollection::TransitiveGroupCollection(Vertex* pVertex, 
                                                     EdgeDir dir) : m_pVertex(pVertex), m_dir(dir)
{

}

//
TransitiveGroup& TransitiveGroupCollection::createGroup(Edge* pIrreducible)
{
	m_groups.push_back(TransitiveGroup(m_pVertex, pIrreducible));
	return m_groups.back();
}

//
TransitiveGroup& TransitiveGroupCollection::getGroup(size_t idx)
{
	assert(idx < m_groups.size());
	return m_groups[idx];
}

// Find the group that the edge with EdgeDesc is in
size_t TransitiveGroupCollection::findGroup(EdgeDesc ed) const
{
	for(size_t i = 0; i < numGroups(); ++i)
	{
		if(m_groups[i].getIrreducible()->getDesc() == ed || m_groups[i].hasTransitive(ed))
			return i;
	}
	assert(false);
	return 0;
}

//
size_t TransitiveGroupCollection::numGroups() const
{
	return m_groups.size();
}

// 
void TransitiveGroupCollection::print() const
{
	for(size_t i = 0; i < numGroups(); ++i)
	{
		std::cout << "Group: " << i << "\n";
		m_groups[i].printMultiOverlap();
	}
}

