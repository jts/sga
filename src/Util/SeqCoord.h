//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SeqCoord - A data structure holding the coordinate
// of a substring of a sequence which consists of an interval and
// the length of the string. Used to build matches and overlaps.
//
#ifndef SEQCOORD_H
#define SEQCOORD_H

#include "Util.h"
#include "Interval.h"

struct SeqCoord
{
	// constructor
	SeqCoord() {}
	SeqCoord(int s, int e, int l) : interval(s, e), seqlen(l) { assert(s <= e); }

	// functions
	inline bool isLeftExtreme() const
	{
		return (interval.start == 0);
	}

	inline bool isRightExtreme() const
	{
		return (interval.end + 1 == seqlen);
	}

	inline bool isExtreme() const
	{
		return (isLeftExtreme() || isRightExtreme());
	}
	
	inline bool isContained() const
	{
		return (isLeftExtreme() && isRightExtreme());
	}

	// Return the length of the interval, which is inclusive
	inline int length() const 
	{ 
		return interval.end - interval.start + 1;
	}

	// Return the distance from the start to the left end
	inline int getLeftDist() const
	{
		return interval.start;
	}

	// Return the distance to the right end
	inline int getRightDist() const
	{
		return seqlen - interval.end - 1;
	}

	// Ensure the interval is within in valid range [0, seqlen)
	inline bool isValid() const
	{
		return interval.start >= 0 && interval.end < seqlen;
	}

	// Flip mirrors the coordinates so they are on the other strand
	// The coordinates are naturally reversed to indicate its the other strand
	inline void flip()
	{
		int tmp = seqlen - interval.start - 1;
		interval.start = seqlen - interval.end - 1;
		interval.end = tmp;
		assert(interval.start <= interval.end);
	}

	// Flip a single position p to the reverse strand for a sequence of length l
	static inline int flip(int p, int l)
	{
		return l - p - 1;
	}

	SeqCoord complement() const;

	// Get the substring described by the interval
	std::string getSubstring(const std::string& str) const;

	// Get the substring described by the complement of the interval
	std::string getComplementString(const std::string& str) const;

	friend std::ostream& operator<<(std::ostream& out, const SeqCoord& sc);
	friend std::istream& operator>>(std::istream& in, SeqCoord& sc);	

	// data
	Interval interval;
	int seqlen;
};

#endif
