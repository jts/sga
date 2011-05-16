//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SampledSuffixArray - Data structure holding
// a subset of suffix array positions. From the
// sampled positions, the other entries of the
// suffix array can be calculated.
//
#ifndef SAMPLED_SUFFIX_ARRAY
#define SAMPLED_SUFFIX_ARRAY

#include "SuffixArray.h"
#include "BWT.h"

class SampledSuffixArray
{
    public:

        SampledSuffixArray();
        SampledSuffixArray(const std::string& filename);
        
        // Calculate the suffix array element for the given index
        SAElem calcSA(int64_t idx, const BWT* pBWT) const;

        // Construct the sampled SA using the initial suffix array index
        // and the bwt of a set of reads
        void build(std::string bwtFilename, std::string saiFilename, int sampleRate = DEFAULT_SA_SAMPLE_RATE);

        // Validate using the full suffix array for the given set of reads. Very slow.
        void validate(std::string readsFile);


    private:

        // SAElems indicating the start of every read in the
        // sequence collection. These elements are in lexicographic order
        // based on the whole read sequence. Tracing a read backwards through
        // the suffix array necessarily ends at one of these positions.
        SAElemVector m_saLexoIndex;

        static const int DEFAULT_SA_SAMPLE_RATE = 64;
        int m_sampleRate;
        SAElemVector m_saSamples;
};

#endif
