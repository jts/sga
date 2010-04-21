//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapAlgorithm - This class implements all the logic
// for finding and outputting overlaps for sequence reads
//
#ifndef OVERLAPALGORITHM_H
#define OVERLAPALGORITHM_H

#include "BWT.h"
#include "OverlapBlock.h"
#include "OverlapSeed.h"
#include "BWTAlgorithms.h"
#include "Util.h"

struct OverlapResult
{
	OverlapResult() : isSubstring(false) {}
	bool isSubstring;
};

class OverlapAlgorithm
{
	public:
		OverlapAlgorithm(const BWT* pBWT, const BWT* pRevBWT, 
		                 size_t minO, double er, bool irrOnly) : m_pBWT(pBWT), 
	                                                            m_pRevBWT(pRevBWT),
																m_minOverlap(minO),
																m_errorRate(er),
																m_bIrreducible(irrOnly) {}

		// Perform the overlap
		// This function is threaded so everything must be const
		OverlapResult overlapRead(const SeqRecord& read, OverlapBlockList* pOutList) const;
	
		// Perform an irreducible overlap
		OverlapResult overlapReadIrreducible(const SeqRecord& read, OverlapBlockList* pOBOut) const;

		// Perform an exhaustive overlap
		OverlapResult overlapReadExhaustive(const SeqRecord& read, OverlapBlockList* pOBOut) const;

		// Perform an inexact overlap
		OverlapResult overlapReadInexact(const SeqRecord& read, OverlapBlockList* pOBOut) const;

		// Write the result of an overlap to an ASQG file
		void writeResultASQG(std::ostream& writer, const SeqRecord& read, const OverlapResult& result) const;

		// Write all the overlap blocks pList to the filehandle
		void writeOverlapBlocks(std::ostream& writer, size_t readIdx, const OverlapBlockList* pList) const;

	private:

		// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
		// overlaps with a suffix of w.
		static void findOverlapBlocks(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
		                              int minOverlap, const AlignFlags& af, OverlapBlockList* pOBTemp, 
							          OverlapBlockList* pOBFinal, OverlapResult& result);


		// Perform an inexact suffix overlap, allowing at most error_rate difference
		static int overlapSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                        double error_rate, int minOverlap, const AlignFlags& af, 
										OverlapBlockList* pOBList, OverlapBlockList* pOBFinal,
										OverlapResult& result);

		static int _alignSegmentSimple(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                       double error_rate, int minOverlap, const AlignFlags& af, 
									   OverlapBlockList* pOBList, OverlapBlockList* pOBFinal, 
									   OverlapResult& result);

		// Using the vector of OverlapBlocks, calculate the irreducible hits
		static void calculateIrreducibleHits(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList* pOBList, OverlapBlockList* pOBFinal);
		
		// Extend each block in obl until all the irreducible overlaps have been found. 
		static void processIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, OverlapBlockList* pOBFinal);

		// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
		static void updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, char b);

		//
        static inline int createOverlapSeeds(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                      double error_rate, int minOverlap, OverlapSeedList* pOutList);

		static inline void extendSeedsExactRight(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                                 ExtendDirection dir, const OverlapSeedList* pInList, 
												 OverlapSeedList* pOutList);

		static inline void extendSeedInexactRight(OverlapSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
		                                          OverlapSeedList* pOutList);

        static inline void extendSeedInexactLeft(OverlapSeed& seed, const std::string& w, const BWT* pBWT, const BWT* pRevBWT,
                                                 OverlapSeedList* pOutList);

		// Data
		const BWT* m_pBWT;
		const BWT* m_pRevBWT;
		size_t m_minOverlap;
		double m_errorRate;
		bool m_bIrreducible;
		

};

#endif
