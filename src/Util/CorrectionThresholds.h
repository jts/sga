//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// CorrectionThresholds.h - Cutoff values for determining
// whether a particular base should be corrected or not
//
#ifndef CORRECTION_THRESHOLDS_H
#define CORRECTION_THRESHOLDS_H

namespace CorrectionThresholds
{
    // The number of reads that needs to support
    // a base call given a low or high quality value
    const static int minSupportLowQuality = 4;
    const static int minSupportHighQuality = 3;

    // The threshold for a determining whether
    // a quality score is low or high
    const static int highQualityCutoff = 20;

};

#endif
