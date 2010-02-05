//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// IntervalSet - Non-overlapping collection of intervals
//
#include "IntervalSet.h"

IntervalSet::IntervalSet()
{

}

// Pre/post-condition: All blocks in the set are non-overlapping and sorted by the lower coordinate (and therefore also the right coordinate)
void IntervalSet::addBlock(const OverlapBlock& block)
{
	(void)block;
	
	// Find all the blocks overlapping this block
	
}
