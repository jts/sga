//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// STCommon.h - Base classes and data structures
//
#ifndef STCOMMON_H
#define STCOMMON_H
#include "STGlobals.h"
//
// Functions
//

// Convert a base to an index
AIdx base2Idx(char b);

// Print out a map using cout
template<class K, class V>
void printMap(const std::map<K,V>& m);

// Print a vector
template<class T>
void printVector(const std::vector<T>& v);

//
// Classes
//

//
// A GSuffix (generalized suffix) is an ID giving the string this suffix belongs to and a number indicating
// the position of the suffix in the string
//
class GSuffix
{
	public:

		// Constructors
		GSuffix() : m_label(0), m_idx(0) {}
		GSuffix(Label l, int i) : m_label(l), m_idx(i) {}
		GSuffix(Label l) : m_label(l), m_idx(0) {}
		
		//
		Label getLabel() const { return m_label; }
		uint16_t getIdx() const { return m_idx; }

		//
		size_t getByteSize() const { return sizeof(m_idx) + sizeof(m_label); }

		// Output
		friend std::ostream& operator<<(std::ostream& out, const GSuffix& gs);

	private:
		Label m_label;
		uint16_t m_idx;
};

//
// A suffix string is a label and the rotated string that represents it
//
class SuffixString
{	
	public:
	
		// Constructors
		SuffixString(Label l, int i, std::string s) : id(l,i), str(s) {}
		SuffixString(Label l, std::string s) : id(l), str(s) {}
		
		// Comparator
		friend int operator<(const SuffixString& o1, const SuffixString& o2);
		
		// Output
		friend std::ostream& operator<<(std::ostream& out, const SuffixString& s);

		// These fields are intentially public
		GSuffix id;
		std::string str;
};

//
// A simple class holding the count for each base of a DNA string (plus the terminator)  
//
typedef uint32_t BaseCount;
class AlphaCount
{
	public:
		AlphaCount();
		void set(char b, BaseCount v);
		void increment(char b);
		BaseCount get(char b) const;
		
		friend std::ostream& operator<<(std::ostream& out, const AlphaCount& ac);

	private:
		BaseCount m_counts[ALPHABET_SIZE];
};

// Typedefs of STL collections of the above classes
typedef std::vector<SuffixString> SuffixStringVector;
typedef std::vector<GSuffix> GSuffixVector;

#endif
