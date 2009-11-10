//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// OverlapQuality - Probabilitic inference ofthe quality of a set of overlaps
//
#ifndef OVERLAPQUALITY_H
#define OVERLAPQUALITY_H

#include "Util.h"
#include "ReadTable.h"
#include "Match.h"
#include "STGlobals.h"

typedef std::vector<unsigned int> uint_vec;
typedef std::vector<double> double_vec;
typedef std::vector<double_vec> double_matrix;
typedef std::vector<std::string> StringVector;

class OverlapQuality
{

	public:
		OverlapQuality(const SeqItem& base_read, const OverlapVector& overlap, const ReadTable* pRT);
		double cluster();

	private:

		// functions

		// Calculate the likelihood of the alignment between x(true sequence) and c
		double calcAlignL(const std::string& x, const std::string& c, const double_vec& r) const;

		// Calculate the expected probability of the match for c
		double calcTrueL(const std::string& c, const double_vec& r) const;

		// Calculate the likelihood of the data given that the actual sequence is x
		double calcLikelihood(const::std::string& x) const;

		// Calculate the most likely base at each position of the data using all the strings
		// that that have index k
		std::string posteriorDecode(size_t k, uint_vec& idx_vec);

		// Calculate the likelihood of the data given the clusterting
		double calculateClusterLikelihood(const StringVector& t_vec, const uint_vec& idx_vec);
	
		// Mask out all the positions of str that match m_baseSeq
		std::string getMaskedSeq(const std::string& str) const;

		// data

		// the sequence of the base read, the read all others were aligned to
		std::string m_baseSeq; 

		// the quality values for the base read 
		double_vec m_baseQual; 		

		// the sequences of all the reads aligned to base, transformed to base's coordinates
		StringVector m_cSeqs; 

		// the quality values for all the reads aligned to base
		double_matrix m_cQuals;
};


#endif
