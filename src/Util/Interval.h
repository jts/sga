//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Interval - A pair of integers denoting a closed interval
//

#ifndef INTERVAL_H
#define INTERVAL_H

#include "Util.h"

struct Interval
{
	// constructors
	Interval() : start(0), end(0) {}
	Interval(int s, int e) : start(s), end(e) {}

	// functions
	friend std::ostream& operator<<(std::ostream& out, const Interval& i);
	friend std::istream& operator>>(std::istream& in, Interval& i);	
	
	// data members
	int start;
	int end;
};

#endif
