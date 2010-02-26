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
#include "Pileup.h"
#include "DNADouble.h"

class MultiOverlap
{
	struct MOData
	{
		MOData(const std::string& s, const Overlap& o) : seq(s), ovr(o), offset(0), partitionID(0), score(0.0f) {}
		static bool sortOffset(const MOData& a, const MOData& b);
		static bool sortID(const MOData& a, const MOData& b);
		
		// data
		std::string seq;
		Overlap ovr;
		int offset;
		int partitionID;
		double score;

		static bool compareScore(const MOData& a, const MOData& b)
		{
			return a.score > b.score;
		}
	};

	typedef std::vector<MOData> MODVector;

	public:

		MultiOverlap(const std::string& rootID, const std::string& rootSeq);
		
		void add(const std::string& seq, const Overlap& ovr);
		void add(const MOData& mod);

		Overlap getOverlap(size_t idx) const;
		size_t getNumBases() const;
		int getPartition(size_t idx) const;
		void setPartition(size_t idx, int p);
		
		// Partition the multioverlap into groups
		void partitionMP(double p_error);
		void partitionLI(double p_error);
		void partitionSL(double p_error);
		void partitionBest(double p_error, size_t n);

		// Count the number of members in the given partition
		size_t countPartition(int id) const;

		double calculateLikelihood() const;
		double calculateGroupedLikelihood() const;
		std::string calculateConsensusFromPartition(double p_error);

		DNADouble calcAlphaProb(size_t idx) const;
		AlphaCount calcAlphaCount(size_t idx) const;
		void calcProb() const;

		// IO
		void print(int default_padding = DEFAULT_PADDING, int max_overhang = DEFAULT_MAX_OVERHANG);
		void printPileup();
		void printGroups();

	private:

		
		Pileup getPileup(int idx) const;
		Pileup getPileup(int idx, int numElems) const;
		Pileup getSingletonPileup(int base_idx, int ovr_idx) const;

		void getPartitionedPileup(int idx, Pileup& g0, Pileup& g1) const;
		PileupVector getPartitionedPileup(int idx, int num_parts) const;

		void printRow(int default_padding, int max_overhang, int root_len, 
		              int offset, int overlap_len, int pid, double score,
					  const std::string& seq, const std::string& id);

		// data
		static const int DEFAULT_PADDING = 20;
		static const int DEFAULT_MAX_OVERHANG = 3;

		std::string m_rootID;
		std::string m_rootSeq;

		MODVector m_overlaps;
};

#endif
