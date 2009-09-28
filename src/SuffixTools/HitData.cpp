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
	uint64_t readID = h.saElem.getID();

	// Check if a hit to this read already exists
	HitMap::iterator iter = m_hitMap.find(readID);

	// If there are no hits, create the vector and add this hit
	if(iter == m_hitMap.end())
	{
		m_hitMap[readID].push_back(h);
	}
	else
	{
		// There are hits to this read stored
		// If the new hit is better, clear out the old hits and append
		// All the hits in the vector have the same score so we can just check the first one
		HitVector& currHits = iter->second;
		assert(currHits.size() > 0);
		if(h.len > currHits.front().len)
		{
			currHits.clear();
			currHits.push_back(h);
		}
		else if(h.len == currHits.front().len)
		{
			currHits.push_back(h);
		}
		// the hit is worse than the hits already stored, do nothing
	}
}

// Turn the hits into a vector
HitVector HitData::getHits() const
{
	HitVector hv;
	hv.reserve(m_hitMap.size());
	for(HitMap::const_iterator iter = m_hitMap.begin(); iter != m_hitMap.end(); ++iter)
	{
		hv.insert(hv.end(), iter->second.begin(), iter->second.end());
	}
	return hv;
}
