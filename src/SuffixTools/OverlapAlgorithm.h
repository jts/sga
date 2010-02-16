//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapAlgorithm - Wrapper class used to find the overlaps
// for a given read using the BWT/FM-index
//
#ifndef OVERLAPALGORITHM_H
#define OVERLAPALGORITHM_H

#include "BWT.h"
#include "OverlapBlock.h"
#include "OverlapSeed.h"
#include "BWTAlgorithms.h"
#include "Util.h"

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
		void overlapRead(const SeqItem& read, OverlapBlockList* pOutList) const;
	
		// Perform an irreducible overlap
		void overlapReadIrreducible(const SeqItem& read, OverlapBlockList* pOBOut) const;

		// Perform an exhaustive overlap
		void overlapReadExhaustive(const SeqItem& read, OverlapBlockList* pOBOut) const;

		// Perform an inexact overlap
		void overlapReadInexact(const SeqItem& read, OverlapBlockList* pOBOut) const;

		// Write all the overlap blocks pList to the filehandle
		void writeOverlapBlocks(size_t readIdx, const OverlapBlockList* pList, std::ofstream& writer) const;

	private:

		// Calculate the ranges in pBWT that contain a prefix of at least minOverlap basepairs that
		// overlaps with a suffix of w.
		static void findOverlapBlocks(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
		                              int minOverlap, const AlignFlags& af, OverlapBlockList* pOBTemp, 
							          OverlapBlockList* pOBFinal);


		// Perform an inexact suffix overlap, allowing at most error_rate difference
		static int overlapSuffixInexact(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                        double error_rate, int minOverlap, const AlignFlags& af, 
										OverlapBlockList* pOBList, OverlapBlockList* pOBFinal);


		// Align a subrange of a string against and fm-index while allowing maxDiff errors. Seeded.
		static int _alignSegment(const std::string& w, int block_start, int block_end,
                                 const BWT* pBWT, const BWT* pRevBWT, int maxDiff, 
								 const AlignFlags& af, OverlapBlockList* pOBList, 
								 OverlapBlockList* pOBFinal);

		static int _alignSegmentSimple(const std::string& w, int block_start, int block_end,
                                       const BWT* pBWT, const BWT* pRevBWT, int maxDiff, 
									   const AlignFlags& af, OverlapBlockList* pOBList, 
									   OverlapBlockList* pOBFinal);



		// Using the vector of OverlapBlocks, calculate the irreducible hits
		static void calculateIrreducibleHits(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList* pOBList, OverlapBlockList* pOBFinal);
		
		// Extend each block in obl until all the irreducible overlaps have been found. 
		static void processIrreducibleBlocks(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, OverlapBlockList* pOBFinal);

		// Update the overlap block list with a righthand extension to b, removing ranges that become invalid
		static void updateOverlapBlockRangesRight(const BWT* pBWT, const BWT* pRevBWT, OverlapBlockList& obList, char b);


		static inline void createOverlapSeeds(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                              int block_start, int block_end, int maxDiff, 
											  OverlapSeedList* pOutList);

		static inline void extendSeedsExactRight(const std::string& w, const BWT* pBWT, const BWT* pRevBWT, 
                                                 ExtendDirection dir, const OverlapSeedList* pInList, 
												 OverlapSeedList* pOutList);

		static inline void extendSeedInexactRight(OverlapSeed& seed, const std::string& w, const BWT* /*pBWT*/, 
                                                  const BWT* pRevBWT, ExtendDirection /*dir*/, OverlapSeedList* pOutList);

		static inline void extendSeedInexactLeft(OverlapSeed& seed, const std::string& w, const BWT* pBWT, 
                                                 const BWT* /*pRevBWT*/, ExtendDirection /*dir*/,
                                                 int block_start, int block_end, int maxDiff, const AlignFlags& af,
												 OverlapSeedList* pOutList, OverlapBlockList* pOBPartialList, 
												 OverlapBlockList* pOBFullList);

		// Data
		const BWT* m_pBWT;
		const BWT* m_pRevBWT;
		size_t m_minOverlap;
		double m_errorRate;
		bool m_bIrreducible;
		

};

#endif
