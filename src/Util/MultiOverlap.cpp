//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MultiOverlap.h - Data structure containing a set
// of overlaps for a given read
//
#include <algorithm>
#include <iostream>
#include "MultiOverlap.h"
#include "Alphabet.h"

MultiOverlap::MultiOverlap(const std::string& rootID, const std::string& rootSeq) : m_rootID(rootID), m_rootSeq(rootSeq)
{
	
}

void MultiOverlap::add(const std::string& seq, const Overlap& ovr)
{
	MOData mod(seq, ovr);

	// Swap root read into first position if necessary
	if(mod.ovr.id[0] != m_rootID)
		mod.ovr.swap();
	assert(mod.ovr.id[0] == m_rootID);

	// RC the sequence if it is different orientation than the root
	if(mod.ovr.match.isRC())
	{
		mod.seq = reverseComplement(mod.seq);
		mod.ovr.match.canonize();
	}

	// Initialize the offset value, the amount that a coordinate 
	// for the non-root sequence must be shifted so that
	// the sequences are aligned
	mod.offset = mod.ovr.match.inverseTranslate(0);
	m_overlaps.push_back(mod);
}

//
Overlap MultiOverlap::getOverlap(size_t idx) const
{
	assert(idx < m_overlaps.size());
	return m_overlaps[idx].ovr;
}

// Compute the likelihood of the multioverlap
double MultiOverlap::calculateLikelihood() const
{
	WARN_ONCE("Compute likelihood using fixed error rate");

	double likelihood = 0.0f;
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup pileup = getPileup(i);
		DNADouble ap = pileup.calculateLikelihoodNoQuality(0.01);
		likelihood += ap.marginalize();
	}
	return likelihood;
}

// Calculate the likelihood of the multioverlap given the groups defined by IntVector
// There are only two possible groups, 0 and 1. The root_sequence (index 0 in the vector) is assumed to belong
// to group 1
double MultiOverlap::calculateGroupedLikelihood(const IntVec& groups) const
{
	double likelihood = 0.0f;
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup g0;
		Pileup g1;
		getGroupedPileup(i, groups, g0, g1);
		if(g0.getDepth() > 0)
			likelihood += g0.calculateLikelihoodNoQuality(0.01).marginalize();
		if(g1.getDepth() > 0)
			likelihood += g1.calculateLikelihoodNoQuality(0.01).marginalize();
	}
	return likelihood;
}


// Calculate the probability of the 4 bases given the multi-overlap
// for a single position
DNADouble MultiOverlap::calcAlphaProb(size_t idx) const
{
	Pileup pileup = getPileup(idx);
	return pileup.calculateSimpleAlphaProb();
}

// Get the number of times each base appears at position
// idx in the multi-align
AlphaCount MultiOverlap::calcAlphaCount(size_t idx) const
{
	Pileup pileup = getPileup(idx);
	return pileup.getAlphaCount();
}


// Calculate the probability of this multioverlap
void MultiOverlap::calcProb() const
{
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup pileup = getPileup(i);
		if(pileup.getDepth() > 1)
		{
			char cnsBase = pileup.calculateSimpleConsensus();
			DNADouble ap = pileup.calculateSimpleAlphaProb();
			char refBase = pileup.getBase(0);
			
			if(refBase != cnsBase)
			{
				int ref_count = pileup.getCount(pileup.getBase(0));
				int cons_count = pileup.getCount(cnsBase);
				int depth = pileup.getDepth();
				double rp = ap.get(refBase);
				double cp = ap.get(cnsBase);
				printf("CNS\t%d\t%d\t%d\t%lf\t%lf\n", ref_count, cons_count, depth, rp, cp);
			}
		}
	}

	/*
	int numMismatches = 0;
	int numAlignedBases = 0;
	double errorRate = 0.01;

	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup pileup = getPileup(i);
		if(pileup.getDepth() > 1)
		{
			char consensus = pileup.calculateSimpleConsensus();
			
			// Calculate the number of bases in the pileup that do not match the consensus
			for(size_t j = 0; j < pileup.getDepth(); ++j)
			{
				if(pileup.getBase(j) != consensus)
					++numMismatches;
				++numAlignedBases;
			}
		}
	}

	double actualRate = double(numMismatches) / double(numAlignedBases);
	double expectedMismatches = errorRate * double(numAlignedBases);
	printf("MM\t%d\t%lf\t%d\t%lf\t%lf\n", numMismatches, expectedMismatches, numAlignedBases, actualRate, errorRate);
	*/
}

size_t MultiOverlap::getNumBases() const
{
	size_t count = 0;
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup p = getPileup(i);
		count += p.getDepth();
	}
	return count;
}

// Get the "stack" of bases that aligns to
// a single position of the root seq, including
// the root base
Pileup MultiOverlap::getPileup(int idx) const
{
	Pileup p;
	p.add(m_rootSeq[idx]);

	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		const MOData& curr = m_overlaps[i];
		// translate idx into the frame of the current sequence
		int trans_idx = idx - curr.offset;
		if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
		{
			p.add(curr.seq[trans_idx]);
		}
	}
	return p;
}

// Fill in the pileup groups g0 and g1 based on the IntVector
void MultiOverlap::getGroupedPileup(int idx, const IntVec& groups, Pileup& g0, Pileup& g1) const
{
	g1.add(m_rootSeq[idx]);

	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		const MOData& curr = m_overlaps[i];
		// translate idx into the frame of the current sequence
		int trans_idx = idx - curr.offset;
		if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
		{
			if(groups[i] == 0)
				g0.add(curr.seq[trans_idx]);
			else
				g1.add(curr.seq[trans_idx]);
		}
	}
}

// Print the MultiOverlap to stdout
void MultiOverlap::print(int default_padding, int max_overhang)
{
	std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::sortOffset);
	std::cout << "\nDrawing overlaps for read " << m_rootID << "\n";
	int root_len = int(m_rootSeq.size());
	
	// Print the root row at the bottom
	printRow(default_padding, max_overhang, root_len, 0, root_len, m_rootSeq, m_rootID);

	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		const MOData& curr = m_overlaps[i];
		int overlap_len = curr.ovr.match.getMaxOverlapLength();

		printRow(default_padding, max_overhang, root_len, 
		         curr.offset, overlap_len, curr.seq, curr.ovr.id[1].c_str());
	}	
}

// Print the MultiOverlap groups specified by the IntVec to stdout
void MultiOverlap::printGroups(const IntVec& groups)
{
	for(size_t i = 0; i <= 1; ++i)
	{
		MultiOverlap mo_group(m_rootID, m_rootSeq);
		for(size_t j = 0; j < m_overlaps.size(); ++j)
		{
			if(groups[j] == (int)i)
				mo_group.add(m_overlaps[j].seq, m_overlaps[j].ovr);
		}
		std::cout << "MO GROUP " << i << "\n";
		mo_group.print();
	}
}


// Print a single row of a multi-overlap to stdout
void MultiOverlap::printRow(int default_padding, int max_overhang, 
                            int root_len, int offset, int overlap_len, 
							const std::string& seq, const std::string& id)
{
	int c_len = seq.length();

	// This string runs from c_offset to c_offset + len
	// Clip the string at -max_overhang to root_len + max_overhang
	int left_clip = std::max(offset, -max_overhang);
	int right_clip = std::min(offset + c_len, root_len + max_overhang);
	
	// translate the clipping coordinates to the string coords
	int t_left_clip = left_clip - offset;
	int t_right_clip = right_clip - offset;
	
	// Calculate the length of the left padding
	int padding = default_padding + left_clip;
	std::string leader = (t_left_clip > 0) ? "..." : "";
	std::string trailer = (t_right_clip < c_len) ? "..." : ""; 
	std::string clipped = seq.substr(t_left_clip, t_right_clip - t_left_clip);
	padding -= leader.size();

	assert(padding >= 0);
	std::string padding_str(padding, ' ');
	std::string outstr = padding_str + leader + clipped + trailer;
	printf("%s\t%d\t%d\t%lf\tID:%s\n", outstr.c_str(), overlap_len, -1, 0.0f, id.c_str());
}

// Print the MultiOverlap horizontally, in a pileup format
void MultiOverlap::printPileup()
{
	std::cout << "\nDrawing overlap pileup for read " << m_rootID << "\n";
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup p = getPileup(i);
		printf("%zu\t%s\n", i, p.toStr().c_str());
	}
}

//
bool MultiOverlap::MOData::sortOffset(const MOData& a, const MOData& b)
{
	return a.offset < b.offset;
}

// Sort by the ID of the non-root sequence, which is 
// guarenteed to bethe second id in the overlap structure
bool MultiOverlap::MOData::sortID(const MOData& a, const MOData& b)
{
	assert(a.ovr.id[0] == b.ovr.id[0]);
	return a.ovr.id[1] < b.ovr.id[1];
}

