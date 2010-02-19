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
#include "Alphabet.h"

#define DEFAULT_PROB 0.01

struct PUElem
{
	PUElem() : base('A'), lp(0) {}
	PUElem(char b, double l) : base(b), lp(l) {}

	char base;

	// log probability this base is correct
	double lp; 
};

typedef std::vector<PUElem> PUElemVector;

class Pileup
{
	public:
		Pileup() {}
		Pileup(size_t n) { m_data.reserve(n); }

		char calculateSimpleConsensus() const;
		AlphaCount getAlphaCount() const;
		AlphaProb calculateSimpleAlphaProb() const;

		void add(char b);
		void add(char b, double lp);

		char getBase(size_t idx) const;
		char getCount(char base) const;
		size_t getDepth() const;

		std::string toStr() const;

	private:
		PUElemVector m_data;
};

typedef std::vector<Pileup> PileupVector;

#endif
