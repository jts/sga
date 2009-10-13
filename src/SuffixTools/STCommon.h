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
#include <utility>
#include <stdint.h>

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
// A Generalized SuffixArray ID (SAElem) is a single number where the high n bits represents the
// identifier of the string (as the index into a StringDictionary) and the low (64 - n) bits 
// represents the position in that string
//
struct SAElem
{
	public:
		SAElem() : m_val(std::numeric_limits<uint64_t>::max()) { }
		SAElem(uint64_t i);
		SAElem(uint64_t i, uint64_t p);

		//
		inline bool isEmpty() const
		{
			return m_val == std::numeric_limits<uint64_t>::max();
		}

		// set the id
		inline void setID(uint64_t i)
		{
			// Clear the HIGH bits by ANDing with the low mask
			m_val &= LOW_MASK;
			
			// Shift the new position into place and set the new value
			i <<= POS_BITS;
			m_val |= i;
		}

		// set the position
		void setPos(uint64_t i)
		{
			// Clear the LOW bits by anding with the high mask
			m_val &= HIGH_MASK;

			// Set the new value
			m_val |= i;
		}

		inline uint64_t getID() const
		{
			return (m_val & HIGH_MASK) >> POS_BITS;
		}

		inline uint64_t getPos() const
		{
			return (m_val & LOW_MASK);
		}

		// Returns true if the suffix is the full length of the string
		inline bool isFull() const
		{
			return getPos() == 0;
		}

		// Input/Output
		friend std::istream& operator>>(std::istream& in, SAElem& s);
		friend std::ostream& operator<<(std::ostream& out, const SAElem& s);


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
		SAElem id;
		std::string str;
};

//
// A simple class holding the count for each base of a DNA string (plus the terminator)  
//
typedef uint64_t BaseCount;
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
typedef std::vector<SAElem> SAElemVector;
typedef std::pair<SAElem, SAElem> SAElemPair;
typedef std::vector<SAElemPair> SAElemPairVec;
typedef std::set<uint64_t> NumericIDSet;

#endif
