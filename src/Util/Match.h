//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Match - A pair of coordinates representing the overlapping
// regions of two sequences
//
#ifndef MATCH_H
#define MATCH_H

#include "Util.h"

struct Match
{
	// functions
	Match() {}
	Match(const SeqCoord& sc1, const SeqCoord& sc2, bool isRC, int nd);
	Match(int s1, int e1, int l1, int s2, int e2, int l2, bool isRC, int nd);

	// Accessors
	inline bool isRC() const { return isReverse; }
	inline int getNumDiff() const { return numDiff; }

	// Translate the SeqCoord c from the frame of coord[0] to coord[1]
	SeqCoord translate(const SeqCoord& c) const;

	// Translate the SeqCoord c from the frame of coord[1] to coord[0]
	SeqCoord inverseTranslate(const SeqCoord& c) const;

	// Translate a single position from c[0] frame to c[1]
	int translate(int c) const;

	// Translate a single position from c[1] frame to c[0]
	int inverseTranslate(int c) const;

	// Return a new match with the coords swapped
	Match swapCoords() const;

	// Swap the coords of this element
	void swap();

	// Flip coord[1] if isReverse is true, effectively
	// bringing the matching strings into the same coordinate system
	void canonize();

	// IO
	friend std::ostream& operator<<(std::ostream& out, const Match& m);
	friend std::istream& operator>>(std::istream& in, Match& m);

	// data
	SeqCoord coord[2];
	bool isReverse;
	int numDiff;
};

// Overlap
struct Overlap
{
	// constructors
	Overlap() {}
	Overlap(std::string i1, const SeqCoord& sc1, std::string i2, const SeqCoord& sc2, bool isRC, int nd); 
	Overlap(std::string i1, int s1, int e1, int l1, std::string i2, int s2, int e2, int l2, bool isRC, int nd); 

	// Swap the order of the elements
	void swap();

	// functions
	friend std::ostream& operator<<(std::ostream& out, const Overlap& o);
	friend std::istream& operator>>(std::istream& in, Overlap& o);

	// data
	std::string id[2];
	Match match;
};

typedef std::vector<Overlap> OverlapVector;

#endif
