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
#include <iostream>

struct Hit
{
	Hit() {}
	Hit(size_t ri, size_t si, uint32_t qs, uint32_t l, bool tr, bool qr) : readIdx(ri), saIdx(si), qstart(qs), 
																			 len(l), targetRev(tr), queryRev(qr) {}

	size_t readIdx; // The index in the read table
	size_t saIdx; // The index into the suffix array

	uint32_t qstart; // The query start position
	uint32_t len; // The overlap length
	bool targetRev; // Whether the target was reversed
	bool queryRev; // Whether the query was reversed

	friend std::ostream& operator<<(std::ostream& out, const Hit& hit)
	{
		out << hit.readIdx << " " << hit.saIdx << " " << hit.qstart << " " << hit.len << " " << hit.targetRev << " " << hit.queryRev;
		return out;
	}

	friend std::istream& operator>>(std::istream& in, Hit& hit)
	{
		in >> hit.readIdx >> hit.saIdx >> hit.qstart >> hit.len >> hit.targetRev >> hit.queryRev;
		return in;
	}

};

typedef std::vector<Hit> HitVector;

#endif
