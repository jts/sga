///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
//
// DindelUtil - Wrappers and utility functions
// for the dindel haplotype scoring functions
//
#include "GraphCompare.h"
#include "DindelRealignWindow.h"

namespace DindelUtil
{
    // Run dindel on a pair of samples
    void runDindelPair(const std::string& normalString, 
                       const std::string& variantString, 
                       const GraphCompareParameters& parameters,
                       VCFFile& baseVCFFile,
                       VCFFile& variantVCFFile);

};
