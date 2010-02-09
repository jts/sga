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

	private:

		// Data
		const BWT* m_pBWT;
		const BWT* m_pRevBWT;
		size_t m_minOverlap;
		double m_errorRate;
		bool m_bIrreducible;
		

};

#endif
