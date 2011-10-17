///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HapgenUtil - Utility functions used for 
// working with the haplotype generation
// modules and output
//
#ifndef HAPGENUTIL_H
#define HAPGENUTIL_H

#include "BWT.h"
#include "SampledSuffixArray.h"

// Simple alignment object representing
// the placement of a string onto a reference genome
struct HapgenAlignment
{

    //
    HapgenAlignment(const int& id, int p, int l, bool rc) : referenceID(id), position(p), length(l), isRC(rc) {}

    //
    int referenceID;
    int position;
    int length; // alignment length
    bool isRC;
};
typedef std::vector<HapgenAlignment> HapgenAlignmentVector;

//
namespace HapgenUtil
{


    // Align the haplotype to the reference genome represented by the BWT/SSA pair
    std::vector<HapgenAlignment> alignHaplotypeToReference(const std::string& haplotype,
                                                           const BWT* pReferenceBWT,
                                                           const SampledSuffixArray* pReferenceSSA);

    // Print an alignment to a reference
    void printAlignment(const std::string& query, const HapgenAlignment& aln, const ReadTable* pRefTable);

};

#endif
