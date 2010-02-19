//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiOverlap.h - Data structure containing a set
// of overlaps for a given read
//
#ifndef MULTIOVERLAP_H
#define MULTIOVERLAP_H

#include "Match.h"

class MultiOverlap
{
	struct MOData
	{
		MOData(const std::string& s, const Overlap& o) : seq(s), ovr(o), offset(0) {}
		static bool sortOffset(const MOData& a, const MOData& b);
		static bool sortID(const MOData& a, const MOData& b);
		
		// data
		std::string seq;
		Overlap ovr;
		int offset;
	};

	typedef std::vector<MOData> MODVector;

	public:

		MultiOverlap(const std::string& rootID, const std::string& rootSeq);
		
		void add(const std::string& seq, const Overlap& ovr);
		void calcProb() const;

		// IO
		void print(int default_padding = DEFAULT_PADDING, int max_overhang = DEFAULT_MAX_OVERHANG);
		void printPileup();

	private:

		
		std::string getPileup(int idx) const;
		char calculateConsensus(const std::string& pileup) const;

		void printRow(int default_padding, int max_overhang, int root_len, 
		              int offset, int overlap_len, const std::string& seq, 
					  const std::string& id);

		// data
		static const int DEFAULT_PADDING = 20;
		static const int DEFAULT_MAX_OVERHANG = 3;

		std::string m_rootID;
		std::string m_rootSeq;

		MODVector m_overlaps;
};

#endif
