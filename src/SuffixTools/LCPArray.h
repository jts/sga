//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// LCPArray - Longest Common Prefix array for a suffix array
// Entry i in the array is the length of the LCP
// between SA entry i and i + 1
//
#ifndef LCPARRAY_H
#define LCPARRAY_H
#include "SuffixArray.h"

typedef std::vector<unsigned int> UIntVector;

class LCPArray
{
    public:
        
        // constructors
        LCPArray(const SuffixArray* pSA, const ReadTable* pRT);

        // getters
        unsigned int get(size_t idx) const;

        //
        size_t countPrefixLength(std::string s1, std::string s2) const;

        // io
        void print(SuffixArray* pSA, ReadTable* pRT) const;

    private:
        UIntVector m_data;
};

#endif
