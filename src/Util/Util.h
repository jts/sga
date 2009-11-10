//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Util - Common data structures and functions
//
#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <istream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include "DNAString.h"

#define CAF_SEP ':'
#define FUNCTION_TIMER Timer functionTimer(__PRETTY_FUNCTION__);
#define VALIDATION_WARNING(x) static bool validation_warn = true; if(validation_warn) \
							   printf("[%s] Warning validation is on\n", (x)); validation_warn = false;

//
// Enums
//
enum EdgeDir
{
        ED_SENSE = 0,
        ED_ANTISENSE,
        ED_COUNT
};

enum EdgeComp
{
        EC_SAME = 0,
        EC_REVERSE
};

const EdgeDir EDGE_DIRECTIONS[ED_COUNT] = { ED_SENSE, ED_ANTISENSE };

//
// Typedef
//
typedef std::vector<int> IntVec;
typedef std::vector<double> DoubleVec;
typedef std::vector<std::string> StringVec;
typedef std::string Sequence;
typedef std::string ContigID;
typedef std::vector<Sequence> SequenceVector;

// SeqItem
struct SeqItem
{
	// data
	std::string id;
	DNAString seq;
};

// KAlignment
struct KAlignment
{
	// functions
	int contigOuterCoordinate() const;
	void flipAlignment(int targetLength); // transform an alignment into the alignment on the rc seq
	int getDistanceToEnd(int targetLen) const; // distance from the align to the end
	static int compareReadPos(const KAlignment& a1, const KAlignment& a2);
	friend std::istream& operator>> (std::istream& in, KAlignment& a);

	// data
	ContigID contig_id;
	int contig_start_pos;
	int read_start_pos;
	int align_length;
	int read_length;
	bool is_reverse;
	
};

// AlignPair
struct AlignPair
{
	// functions
	friend std::istream& operator>> (std::istream& in, AlignPair& ap)
	{
		std::string readname;
		in >> readname >> ap.aligns[0] >> readname >> ap.aligns[1];
		return in;
	}

	// data
	KAlignment aligns[2];
};

// AdjInfo
struct AdjInfo
{

	// functions
	friend std::istream& operator>>(std::istream& in, AdjInfo& a);

	// data
	ContigID from;
	ContigID to;
	int dir;
	bool comp;
};

// Interval 
struct Interval
{
	// constructors
	Interval() : start(0), end(0) {}
	Interval(int s, int e) : start(s), end(e) {}

	// functions
	friend std::ostream& operator<<(std::ostream& out, const Interval& i);
	friend std::istream& operator>>(std::istream& in, Interval& i);	
	
	// data 
	int start;
	int end;
};

// String, coordinate pair
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

//
// Functions
//
std::string stripFilename(std::string filename);
void checkFileHandle(std::ifstream& fh, std::string fn);

// Key-value operations
template <class C>
std::string makeKeyValue(std::string key, C value)
{
	std::stringstream ss;
	ss << key << CAF_SEP << value;
	return ss.str();
}

StringVec split(std::string in, char delimiter);
void splitKeyValue(std::string in, std::string& key, std::string& value);


// Debug function to get the distance between two reads based on their names, which 
// encodes the positions
size_t debug_getReadDistFromNames(const std::string& name1, const std::string& name2);

//
// Return the lexographic value for the given base
//
static const uint8_t s_lexoRankLUT[256] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
	0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

inline static uint8_t getBaseRank(char b)
{
	return s_lexoRankLUT[static_cast<uint8_t>(b)];
}

//
// Sequence operations
//
Sequence reverseComplement(const Sequence& seq);
Sequence complement(const Sequence& seq);
Sequence reverse(const Sequence& seq);

// Count the number of differences between s1 and s2 over the first n chars
int countDifferences(const std::string& s1, const std::string& s2, size_t n);


// Complement a base
inline char complement(char base)
{
	switch(base)
	{
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		default:
			assert(false && "Unknown base!");
			return 'N';
	}
}


//
// Edge Operations
//
inline EdgeDir operator!(const EdgeDir& dir)
{
        return (dir == ED_SENSE) ? ED_ANTISENSE : ED_SENSE;
}

inline EdgeComp operator!(const EdgeComp& comp)
{
        return (comp == EC_SAME) ? EC_REVERSE : EC_SAME;
}

// Correct an edges direction, given the relationship between the nodes
inline EdgeDir correctDir(EdgeDir dir, EdgeComp comp)
{
	return (comp == EC_SAME) ? dir : !dir;
}

#endif

