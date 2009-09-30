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

#define CAF_SEP ':'

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
	std::string id;
	Sequence seq;
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
	SeqCoord(std::string i, int s, int e, int l) : id(i), interval(s, e), seqlen(l) {}

	// functions
	bool isLeftExtreme() const;
	bool isRightExtreme() const;
	bool isExtreme() const;
	bool isContained() const;
	bool isReverse() const;

	// Get the substring described by the interval
	std::string getSubstring(std::string str) const;

	friend std::ostream& operator<<(std::ostream& out, const SeqCoord& sc);
	friend std::istream& operator>>(std::istream& in, SeqCoord& sc);	

	// data
	std::string id;
	Interval interval;
	int seqlen;
};

// Overlap
struct Overlap
{
	// constructors
	Overlap() {}
	Overlap(std::string i1, int s1, int e1, int l1, std::string i2, int s2, int e2, int l2); 

	// functions
	friend std::ostream& operator<<(std::ostream& out, const Overlap& o);
	friend std::istream& operator>>(std::istream& in, Overlap& o);

	// data
	SeqCoord read[2];
};

typedef std::vector<Overlap> OverlapVector;

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

//
// Sequence operations
//
Sequence reverseComplement(const Sequence& seq);
Sequence complement(const Sequence& seq);
Sequence reverse(const Sequence& seq);
char complement(char base);

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

