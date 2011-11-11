//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTCABauerCoxRosone - In-memory version of Illumina's
// BWT construction algorithm
#ifndef BWTCA_COX_BAUER_H

#include "ReadTable.h"



namespace BWTCA
{
    // Element of the N array
    struct NElem
    {
        NElem() : index(0), sym('\0') {}

        // Comparator. Sort first on the symbol, then break ties using the index
        friend bool operator<(const NElem& a, const NElem& b)
        {
            return a.sym < b.sym;
        }

        // Data
        uint32_t index;
        char sym;
    };
    typedef std::vector<NElem> NVector;

    // Construct the burrows-wheeler transform of the table of reads
    void runBauerCoxRosone(const ReadTable* pRT);
    
    // Run the initial special first iteration of the algorithm
    void outputInitialCycle(const ReadTable* pRT, NVector& n_vector, std::string& bwt);
};

#endif
