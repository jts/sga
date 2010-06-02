//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Alphabet.h - Abstraction of the alphabet
// that is used in the suffix array, fm-index,
// etc data structures
//
#ifndef ALPHABET_H
#define ALPHABET_H

#include <utility>
#include <stdint.h>
#include <limits>
#include <math.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <iostream>
#include <algorithm>
#include "Util.h"

//
// Constants
//
// TODO: Refactor these into a namespaced alphabet like DNA_ALPHABET below
const uint8_t ALPHABET_SIZE = 5;
const char ALPHABET[ALPHABET_SIZE] = {'A', 'C', 'G', 'T', '$'};
const char RANK_ALPHABET[ALPHABET_SIZE] = {'$', 'A', 'C', 'G', 'T'};
const uint8_t DNA_ALPHABET_SIZE = 4;
typedef uint64_t BaseCount;

namespace DNA_ALPHABET
{
    static const uint8_t s_dnaLexoRankLUT[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
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
        return s_dnaLexoRankLUT[static_cast<uint8_t>(b)];
    }
    
    inline char getBase(size_t idx)
    {
        assert(idx < DNA_ALPHABET_SIZE);
        return ALPHABET[idx];
    }

    static const uint8_t size = 4;
};

namespace BWT_ALPHABET
{
    static const uint8_t s_bwtLexoRankLUT[256] = {
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

    static const uint8_t size = 5;

    inline static uint8_t getRank(char b)
    {
        return s_bwtLexoRankLUT[static_cast<uint8_t>(b)];
    }
    
    inline char getChar(size_t idx)
    {
        return RANK_ALPHABET[idx];
    }
};

//
// A simple class holding the count for each base of a DNA string (plus the terminator)
// Note that this uses RANK_ALPHABET
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

        //
        inline char getMaxBase() const
        {
            char base;
            BaseCount val;
            getMax(base, val);
            return base;
        }

        //
        inline BaseCount getMaxCount() const
        {
            char base;
            BaseCount val;
            getMax(base, val);
            return val;
        }

        //
        inline void getMax(char& base, BaseCount& val) const
        {
            BaseCount max = 0;
            int maxIdx = 0;
            for(int i = 0; i < ALPHABET_SIZE; ++i)
            {
                if(m_counts[i] > max)
                {
                    max = m_counts[i];
                    maxIdx = i;
                }
            }
            base = RANK_ALPHABET[maxIdx];
            val = max;
        }

        // Sort the DNA bases into count order and
        // write them to pOut. This does not include the '$' symbol
        inline void getSorted(char* pOut, size_t len)
        {
            assert(len >= ALPHABET_SIZE);
            (void)len;
            for(size_t i = 0; i < ALPHABET_SIZE; ++i)
                pOut[i] = RANK_ALPHABET[i];
            std::sort(pOut, pOut+ALPHABET_SIZE, AlphaCountCompareDesc(this));
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

        struct AlphaCountCompareDesc
        {
            AlphaCountCompareDesc(AlphaCount* p) : pAC(p) {}
            bool operator()(char a, char b)
            {
                return pAC->get(a) > pAC->get(b);
            }
            const AlphaCount* pAC;
        };

};

#endif
