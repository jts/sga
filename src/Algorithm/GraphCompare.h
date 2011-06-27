///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphCompare - Compare two (abstractly represented)
// graphs against each other to find strings
// only present in one.
//
// The graphs are abstractly represented as
// an FM-index.
//
#ifndef GRAPH_COMPARE_H
#define GRAPH_COMPARE_H

#include <list>
#include "BWT.h"

// Typedefs


class GraphCompare
{
    public:

        //
        // Functions
        //
        GraphCompare(const BWT* pBaseBWT, const BWT* pBaseRBWT,
                     const BWT* pVariantBWT, const BWT* pVariantRBWT,
                     int kmer);

        ~GraphCompare();
        
        // Perform the comparison
        void run();

    private:

        //
        // Functions
        //

        //
        // Data
        //
        const BWT* m_pBaseBWT; 
        const BWT* m_pBaseRevBWT;
        const BWT* m_pVariantBWT; 
        const BWT* m_pVariantRevBWT;
        int m_kmer;
};

#endif
