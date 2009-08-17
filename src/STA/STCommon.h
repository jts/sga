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
// A Generalized SuffixArray ID (SAID) is a single number where the high n bits represents the
// identifier of the string (as the index into a StringDictionary) and the low (64 - n) bits 
// represents the position in that string
//
struct SAID
{
	public:
		SAID() : m_val(0) {}
		SAID(uint64_t i);
		SAID(uint64_t i, uint64_t p);

		//
		void setID(uint64_t i);
		void setPos(uint64_t i);
		uint64_t getID() const;
		uint64_t getPos() const;

		// Output
		friend std::ostream& operator<<(std::ostream& out, const SAID& s);


	private:
		
		//
		uint64_t m_val;

		// Masks
		static const uint8_t ID_BITS = 36; // Allows up to 68 billion IDs
		static const uint8_t POS_BITS = 64 - ID_BITS;
		static const uint64_t HIGH_MASK = ~0 << POS_BITS;
		static const uint64_t LOW_MASK = ~HIGH_MASK;
};


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
		SuffixString(int i, int p, std::string s) : id(i,p), str(s) {}
		SuffixString(int i, std::string s) : id(i), str(s) {}
		
		// Comparator
		friend int operator<(const SuffixString& o1, const SuffixString& o2);
		
		// Output
		friend std::ostream& operator<<(std::ostream& out, const SuffixString& s);

		// These fields are intentially public
		SAID id;
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

#endif
