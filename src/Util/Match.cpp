//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Match - A pair of coordinates representing the overlapping
// regions of two sequences
//
#include "Match.h"

// Constructor
Match::Match(const SeqCoord& sc1, const SeqCoord& sc2, bool isRC, int nd)
{
	coord[0] = sc1;
	coord[1] = sc2;
	isReverse = isRC;
	numDiff = nd;
}

// Constructor
Match::Match(int s1, int e1, int l1, int s2, int e2, int l2, bool isRC, int nd)
{
	coord[0] = SeqCoord(s1, e1, l1);
	coord[1] = SeqCoord(s2, e2, l2);
	isReverse = isRC;
	numDiff = nd;
}

// Swap the order of the elements
void Match::swap()
{
	SeqCoord tmp = coord[0];
	coord[0] = coord[1];
	coord[1] = tmp;
}

// Flip coord2 so that it is in the same frame as coord1
void Match::canonize()
{
	if(isRC())
		coord[1].flip();
	isReverse = false;
}

// Translate the SeqCoord c from the frame of coord[0] to coord[1]
SeqCoord Match::translate(const SeqCoord& c) const
{
	assert(coord[0].length() == coord[1].length()); // ensure translation is valid
	int t = coord[1].interval.start - coord[0].interval.start;
	
	SeqCoord out;
	out.seqlen = coord[1].seqlen;
	out.interval.start = c.interval.start + t;
	out.interval.end = c.interval.end + t;
	return out;
}

// Translate the SeqCoord c from the frame of coord[1] to coord[0]
SeqCoord Match::inverseTranslate(const SeqCoord& c) const
{
	assert(coord[0].length() == coord[1].length()); // ensure translation is valid
	int t = coord[0].interval.start - coord[1].interval.start;
	
	SeqCoord out;
	out.seqlen = coord[0].seqlen;
	out.interval.start = c.interval.start + t;
	out.interval.end = c.interval.end + t;
	return out;
}

// Translate a single position from c[0] frame to c[1]
int Match::translate(int c) const
{
	assert(coord[0].length() == coord[1].length()); // ensure translation is valid
	int t = coord[1].interval.start - coord[0].interval.start;
	return c + t;
}

// Translate a single position from c[1] frame to c[0]
int Match::inverseTranslate(int c) const
{
	assert(coord[0].length() == coord[1].length());
	int t = coord[0].interval.start - coord[1].interval.start;
	return c + t;
}	

// Output
std::ostream& operator<<(std::ostream& out, const Match& m)
{
	out << m.coord[0] << " " << m.coord[1] << " " << m.isReverse << " " << m.numDiff;
	return out;
}

// Input
std::istream& operator>>(std::istream& in, Match& m)
{
	in >> m.coord[0] >> m.coord[1] >> m.isReverse >> m.numDiff;
	return in;
}

//
// Overlap
//
Overlap::Overlap(std::string i1, const SeqCoord& sc1, std::string i2, const SeqCoord& sc2, bool isRC, int nd)
{
	id[0] = i1;
	id[1] = i2;
	match = Match(sc1, sc2, isRC, nd);
}

Overlap::Overlap(std::string i1, int s1, int e1, int l1,
                 std::string i2, int s2, int e2, int l2, bool isRC, int nd)
{
	id[0] = i1;
	id[1] = i2;
	match = Match(s1, e1, l1, s2, e2, l2, isRC, nd);
}

// Swap the elements
void Overlap::swap()
{
	id[0].swap(id[1]);
	match.swap();
}

// Output
std::ostream& operator<<(std::ostream& out, const Overlap& o)
{
	out << o.id[0] << " " << o.id[1] << " " << o.match;
	return out;
}

// Input
std::istream& operator>>(std::istream& in, Overlap& o)
{
	in >> o.id[0] >> o.id[1] >> o.match;
	return in;
}
