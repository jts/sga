//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// HitData - Data structure holding all positions of
// alignment hits of a read
//
#ifndef OVERLAPDATA_H
#define OVERLAPDATA_H
#include "STCommon.h"

struct Hit
{
	Hit() {}
	Hit(SAID i, uint32_t qs, uint32_t l) : said(i), qstart(qs), len(l) {}
	SAID said;
	uint32_t qstart;
	uint32_t len;
};

typedef std::map<uint64_t, Hit> HitMap;
typedef std::vector<Hit> HitVector;

class HitData
{
	public:
		
		HitData() {}
		void addHit(const Hit& h);
		HitVector getHits() const;

	private:
		HitMap m_hitMap;
};

#endif
