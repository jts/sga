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
#include <stdint.h>
#include "DNAString.h"

#define CAF_SEP ':'
#define FUNCTION_TIMER Timer functionTimer(__PRETTY_FUNCTION__);
#define VALIDATION_WARNING(x) static bool validation_warn = true; if(validation_warn) \
							   printf("[%s] Warning validation is on\n", (x)); validation_warn = false;

#define WARN_ONCE(x) static bool _warn_once = true; if(_warn_once) \
					 printf("WARNING: [%s]\n", (x)); _warn_once = false;


//
// Typedef
//
typedef std::string Sequence;
typedef std::string ContigID;
typedef std::vector<int> IntVec;
typedef std::vector<double> DoubleVec;
typedef std::vector<std::string> StringVec;
typedef std::vector<Sequence> SequenceVector;

// SeqItem
struct SeqItem
{
	// data
	std::string id;
	DNAString seq;
};

//
// Functions
//
std::string stripFilename(std::string filename);
void checkFileHandle(std::ifstream& fh, std::string fn);
void checkFileHandle(std::ofstream& fh, std::string fn);

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
std::string getPairID(const std::string& id);

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
std::string reverseComplement(const std::string& seq);
std::string complement(const std::string& seq);
std::string reverse(const std::string& seq);

// Return the prefix/suffix of seq of length len 
std::string prefix(const std::string& seq, const unsigned int len);
std::string suffix(const std::string& seq, const unsigned int len);


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

#endif
