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
#include <stack>
#include "BWT.h"
#include "BWTInterval.h"

// Typedefs
typedef std::vector<const BWT*> BWTVector;

struct GraphCompareStackNode
{
    // data

    static const size_t NUM_GRAPHS = 2;

    BWTIntervalPair intervalPairs[NUM_GRAPHS];
    AlphaCount64 lowerCounts[NUM_GRAPHS];
    AlphaCount64 upperCounts[NUM_GRAPHS];
    std::string str;
    uint8_t alphaIndex;

    // functions

    // Initialize the intervals to be the range of all strings starting with b
    void initialize(char b, const BWTVector& bwts, const BWTVector& rbwts);

    // Update the intervals to extend to symbol b
    void update(char b, const BWTVector& bwts, const BWTVector& rbwts);
    
    //
    AlphaCount64 getAggregateExtCount() const;

    //
    void print() const;

};
typedef std::stack<GraphCompareStackNode> GraphCompareStack;

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
        size_t m_kmer;
};

#endif
