//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapArray - Data structure and functions used to count 
// the number of times a suffix of a given rank occurs in a data set
//
#ifndef GAPARRAY_H
#define GAPARRAY_H

#include <vector>
#include "BWT.h"
#include "DNAString.h"

typedef uint32_t GAP_TYPE;
typedef std::vector<GAP_TYPE> GapArray;

void updateGapArray(const DNAString& w, const BWT* pBWTInternal, GapArray& gap_array);
void incrementGapArray(int64_t rank, GapArray& gap_array);

#endif
