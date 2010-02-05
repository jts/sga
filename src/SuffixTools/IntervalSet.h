//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// IntervalSet - Non-overlapping collection of intervals
//
#ifndef INTERVALSET_H
#define INTERVALSET_H

#include <list>
#include "OverlapData.h"

typedef std::list<OverlapBlock> BlockSet;

class IntervalSet
{
	public:
		
		IntervalSet();
		void addBlock(const OverlapBlock& block);

	private:

		BlockSet m_set;
};

#endif
