///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringBuilder - Iteratively construct strings
// that represent sequences in an assembly graph.
// The assembly graph is abstractly represented as
// an FM-index.
//
#ifndef STRING_BUILDER_H
#define STRING_BUILDER_H

#include "BWT.h"

class StringBuilder
{
    public:
        StringBuilder(const std::string& seed, int kmer, const BWT* pBWT, const BWT* pRevBWT);

        // Perform one round of extension for all the strings
        void extendOnce();

        // Remove all strings but the top n scoring
        void cull(const std::string& query, size_t n);

        // Print the set of strings created
        void print() const;
        void print(const std::string& query) const;

    private:

        //
        const BWT* m_pBWT; 
        const BWT* m_pRevBWT;

        StringVector m_strings;
        int m_kmer;
};

#endif
