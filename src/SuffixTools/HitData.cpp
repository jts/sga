//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// HitData - Data structure holding all positions of
// alignment hits of a read
//
#include "HitData.h"

// Add a hit to the map
void HitData::addHit(const Hit& h)
{
	uint64_t readID = h.said.getID();

	// Check if a hit to this read already exists
	HitMap::iterator iter = m_hitMap.find(readID);

	// Insert/Replace
	if(iter == m_hitMap.end() || iter->second.len < h.len)
	{
		m_hitMap[readID] = h;
	}
}

// Turn the hits into a vector
HitVector HitData::getHits() const
{
	HitVector hv;
	hv.reserve(m_hitMap.size());
	for(HitMap::const_iterator iter = m_hitMap.begin(); iter != m_hitMap.end(); ++iter)
	{
		hv.push_back(iter->second);
	}
	return hv;
}
