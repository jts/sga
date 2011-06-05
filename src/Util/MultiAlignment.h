//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiAlignment - Class constructing a multiple
// alignment between a root sequence and a set of sequences
//
#ifndef MULTIALIGNMENT_H
#define MULTIALIGNMENT_H

#include <string>
#include <vector>
#include "Util.h"

struct MAlignData
{
    std::string str;
    std::string expandedCigar; // cigar string with symbol counts expanded
    std::string padded;

    int position; // start position on the root
    int targetID; 
    int targetAlignLength;
    int targetPosition;

    // Comparator
    static bool sortPosition(const MAlignData& a, const MAlignData& b) { return a.position < b.position; }

};
typedef std::vector<MAlignData> MAlignDataVector;

class MultiAlignment
{
    public:
        MultiAlignment(std::string rootStr, const MAlignDataVector& inData);

        std::string generateConsensus();
        void print(const std::string* pConsensus = NULL) const;

    private:

        
        // data
        MAlignDataVector m_alignData;
        int m_verbose;
};

#endif
