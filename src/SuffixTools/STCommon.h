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
#include "Util.h"
#include <utility>
#include <stdint.h>
#include <limits>
#include <emmintrin.h>
#include <xmmintrin.h>

//
// Functions
//

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
// A simple class holding the count for each base of a DNA string (plus the terminator)
// Note that the counts are stored in lexographic order
//
typedef uint64_t BaseCount;
class AlphaCount
{
	public:
	
		//
		inline AlphaCount()
		{
			memset(m_counts, 0, ALPHABET_SIZE * sizeof(BaseCount));
		}

		//
		inline void set(char b, BaseCount v)
		{
			m_counts[getBaseRank(b)] = v;
		}

		// 
		inline void increment(char b)
		{
			m_counts[getBaseRank(b)]++;
		}

		// 
		inline BaseCount get(char b) const
		{
			return m_counts[getBaseRank(b)];
		}

		//
		inline BaseCount getByIdx(const int i) const
		{
			return m_counts[i];
		}

		// Return the base for index i
		static char getBase(size_t i)
		{
			return RANK_ALPHABET[i];
		}

		// Swap the (A,T) and (C,G) entries, which turns the AlphaCount
		// into the AlphaCount for the complemented sequence
		inline void complement()
		{
			BaseCount tmp;

			// A,T
			tmp = m_counts[4];
			m_counts[4] = m_counts[1];
			m_counts[1] = tmp;

			// C,G
			tmp = m_counts[3];
			m_counts[3] = m_counts[2];
			m_counts[2] = tmp;
		}

		// Return the sum of the basecounts for characters lexo. lower than b
		inline BaseCount getLessThan(char b) const
		{
			BaseCount out = 0;
			int stop = getBaseRank(b);
			for(int i = 0; i < stop; ++i)
				out += m_counts[i];
			return out;
		}

		// Returns true if only one of the DNA characters
		// has a non-zero count
		inline bool hasUniqueDNAChar()
		{
			// index 0 is the '$' character, which we skip
			bool nonzero = false;
			for(int i = 1; i < ALPHABET_SIZE; ++i)
			{
				if(m_counts[i] > 0)
				{
					// if nonzero is set, there is some other nonzero character, return false
					if(nonzero)
						return false;
					else
						nonzero = true;
				}
			}
			return nonzero;
		}

		// Returns true if any DNA char 
		// has a non-zero count
		inline bool hasDNAChar()
		{
			// index 0 is the '$' character, which we skip
			for(int i = 1; i < ALPHABET_SIZE; ++i)
			{
				if(m_counts[i] > 0)
					return true;
			}
			return false;
		}
		

		// Return the unique DNA character described by the alphacount
		// Returns '$' if no such character exists
		// Asserts if more than one character is described by the alphacount
		inline char getUniqueDNAChar()
		{
			char r = '$';
			for(int i = 1; i < ALPHABET_SIZE; ++i)
			{
				if(m_counts[i] > 0)
				{
					assert(r == '$');
					// lookup the character in the ranked alphabet, since the elements
					// are stored in lexographic order
					r = RANK_ALPHABET[i];
				}
			}
			return r;
		}

		// Operators
		friend std::ostream& operator<<(std::ostream& out, const AlphaCount& ac);
		friend std::istream& operator>>(std::istream& in, AlphaCount& ac);
		
#define USE_SSE 0

		inline friend AlphaCount operator+(const AlphaCount& left, const AlphaCount& right)
		{
			AlphaCount out;

#if USE_SSE	
			__m128i a = _mm_loadu_si128( (__m128i *) &left.m_counts[0]);
			__m128i b = _mm_loadu_si128( (__m128i *) &right.m_counts[0]);
		    __m128i c = _mm_add_epi64(a, b);
			_mm_store_si128( (__m128i *) &out.m_counts[0], c);

			a = _mm_loadu_si128( (__m128i *) &left.m_counts[2]);
			b = _mm_loadu_si128( (__m128i *) &right.m_counts[2]);
		    c = _mm_add_epi64(a, b);
			_mm_store_si128( (__m128i *) &out.m_counts[2], c);
			out.m_counts[4] = left.m_counts[4] + right.m_counts[4];
#else
			out.m_counts[0] = left.m_counts[0] + right.m_counts[0];
			out.m_counts[1] = left.m_counts[1] + right.m_counts[1];
			out.m_counts[2] = left.m_counts[2] + right.m_counts[2];
			out.m_counts[3] = left.m_counts[3] + right.m_counts[3];
			out.m_counts[4] = left.m_counts[4] + right.m_counts[4];
#endif
			return out;
		}

		inline AlphaCount& operator+=(const AlphaCount& other)
		{

#if USE_SSE
			__m128i a = _mm_loadu_si128( (__m128i *) &m_counts[0]);
			__m128i b = _mm_loadu_si128( (__m128i *) &other.m_counts[0]);
		    __m128i c = _mm_add_epi64(a, b);
			_mm_store_si128( (__m128i *) &m_counts[0], c);

			a = _mm_loadu_si128( (__m128i *) &m_counts[2]);
			b = _mm_loadu_si128( (__m128i *) &other.m_counts[2]);
		    c = _mm_add_epi64(a, b);
			_mm_store_si128( (__m128i *) &m_counts[2], c);
#else
			
			m_counts[0] += other.m_counts[0];
			m_counts[1] += other.m_counts[1];
			m_counts[2] += other.m_counts[2];
			m_counts[3] += other.m_counts[3];
			m_counts[4] += other.m_counts[4];
#endif
			return *this;
		}

		// As the counts are unsigned integers, each value in left
		// must be larger or equal to value in right. The calling function
		// must guarentee this.
		friend AlphaCount operator-(const AlphaCount& left, const AlphaCount& right)
		{
			AlphaCount out;
#if USE_SSE
			__m128i a = _mm_loadu_si128( (__m128i *) &left.m_counts[0]);
			__m128i b = _mm_loadu_si128( (__m128i *) &right.m_counts[0]);
		    __m128i c = _mm_sub_epi64(a, b);
			_mm_store_si128( (__m128i *) &out.m_counts[0], c);

			a = _mm_loadu_si128( (__m128i *) &left.m_counts[2]);
			b = _mm_loadu_si128( (__m128i *) &right.m_counts[2]);
		    c = _mm_sub_epi64(a, b);
			_mm_store_si128( (__m128i *) &out.m_counts[2], c);
			out.m_counts[4] = left.m_counts[4] - right.m_counts[4];
#else
			out.m_counts[0] = left.m_counts[0] - right.m_counts[0];
			out.m_counts[1] = left.m_counts[1] - right.m_counts[1];
			out.m_counts[2] = left.m_counts[2] - right.m_counts[2];
			out.m_counts[3] = left.m_counts[3] - right.m_counts[3];
			out.m_counts[4] = left.m_counts[4] - right.m_counts[4];
#endif
			return out;
		}

	private:
		BaseCount m_counts[ALPHABET_SIZE];
};

// Typedefs of STL collections of the above classes
typedef std::vector<SAElem> SAElemVector;
typedef std::pair<SAElem, SAElem> SAElemPair;
typedef std::vector<SAElemPair> SAElemPairVec;
typedef std::set<uint64_t> NumericIDSet;

#endif
