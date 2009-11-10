//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Pileup - An array of strings containing all the seen
// bases for a given region/read
//
#ifndef PILEUP_H
#define PILEUP_H
#include "Util.h"

struct PUElem
{
	PUElem() : base('A'), lp(1) {}
	PUElem(char b, double l) : base(b), lp(l) {}
	char base;
	double lp; // log probability this base is correct
};

typedef std::vector<PUElem> PUElemVector;

class Pileup
{
	public:
		Pileup(size_t n) { m_data.reserve(n); }
		void add(char b, double p);
		std::string getStr() const;

	private:
		PUElemVector m_data;
};

typedef std::vector<Pileup> PileupVector;

#endif
