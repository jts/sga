//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Quality - functions for manipulating quality values
//
#ifndef QUALITY_H
#define QUALITY_H

namespace Quality
{
	inline int char2phred(char b)
	{
		uint8_t v = b;
		assert(v >= 33);
		return v - 33;
	}

	inline char phred2char(int p)
	{
		uint8_t v = (p <= 93) ? p : 93;
		char c = v + 33;
		return c;
	}
};


#endif
