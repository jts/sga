//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapTools - Wrapper for the overlap machinery 
// to perform an overlap computation for two strings
//
#ifndef OVERLAPTOOLS_H
#define OVERLAPTOOLS_H

#include "Match.h"
#include "DPAlignment.h"
#include <list>

namespace OverlapTools
{
    // Datatypes
    typedef std::vector<int> IntVector;
    typedef std::vector<IntVector> DPMatrix;

    //typedef EditDistanceScoring DPScoring;
    typedef SimilarityScoring DPScoring;
    
    void dpOverlap(const std::string& s1, const std::string& s2);

    void findBestOverlap(const std::string& s1, const std::string& s2, 
                         int minOverlap, double maxErrorRate, const DPAlignment& dpAlign);

    void findBestOverlapByScore(const std::string& s1, const std::string& s2, 
                                int minOverlap, double maxErrorRate, const DPAlignment& dpAlign);

};

#endif
