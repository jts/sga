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

//
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
void MultiOverlap::add(const MOData& mod)
{
	assert(mod.ovr.id[0] == m_rootID);
	m_overlaps.push_back(mod);
}

//
Overlap MultiOverlap::getOverlap(size_t idx) const
{
	assert(idx < m_overlaps.size());
	return m_overlaps[idx].ovr;
}

// Returns true if the overlap at idx has the include flag set
int MultiOverlap::getPartition(size_t idx) const
{
	assert(idx < m_overlaps.size());
	return m_overlaps[idx].partitionID;
}

// Returns true if the overlap at idx has the include flag set
void MultiOverlap::setPartition(size_t idx, int p)
{
	assert(idx < m_overlaps.size());
	m_overlaps[idx].partitionID = p;
}

// Return the total number of bases in the multioverlap
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

//
void MultiOverlap::partitionMP(double p_error)
{
	double initial_likelihood = calculateGroupedLikelihood();

	// Compute the likelihood of the alignment between the root sequence
	// and every member of the multioverlap
	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		double likelihood = 0.0f;
		double total_bases = 0;
		for(size_t j = 0; j < m_rootSeq.size(); ++j)
		{
			Pileup pileup = getSingletonPileup(j, i);
			DNADouble ap = pileup.calculateLikelihoodNoQuality(p_error);
			likelihood += ap.marginalize(0.25f);
			total_bases += pileup.getDepth();
		}
		m_overlaps[i].score = likelihood / total_bases;
		
		// Initially set all elements to the "other" partition
		m_overlaps[i].partitionID = 1;
	}

	// Sort the overlaps by score
	std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::compareScore);

	/*
	std::cout << "Print score list\n";
	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		MOData& curr = m_overlaps[i];
		std::cout << i << "\t" << curr.ovr.id[1] << "\t" << curr.score << "\t" << curr.partitionID << "\n";
	}
	*/
	/*
	double prev_likelihood = calculateGroupedLikelihood();
	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		m_overlaps[i].partitionID = 0;
		double likelihood = calculateGroupedLikelihood();
		double improvement = likelihood - prev_likelihood;
		std::cout << "Elem: " << i << " likelihood: " << likelihood << " improvement: " << improvement << "\n";
		prev_likelihood = likelihood;
	}
	*/

	double max = -std::numeric_limits<double>::max();
	size_t num_elems_max = 0;
	size_t num_total = m_overlaps.size();
	double prior = 0.25f;
	double log_prior = log(prior);

	for(size_t j = 0; j < num_total; ++j)
	{
		double post_max_sum = 0.0f;
		for(size_t i = 0; i < m_rootSeq.size(); ++i)
		{
			Pileup pileup = getPileup(i, j);
			DNADouble ap = pileup.calculateLikelihoodNoQuality(p_error);
			double marginal = ap.marginalize(prior);
			//double marginal = ap.marginalize(prior);
			double ratio = log_prior - marginal;
			ap.add(ratio);
			double mv = ap.getMaxVal();
			if(mv > 0.0f)
				mv = 0.0f;
			//std::cout << "postmax: " << ap.getMaxVal() << " clamped: " << mv << "\n";
			post_max_sum += mv;
			//std::cout << "sumnow: " << post_max_sum << "\n";
		}
		
		if(post_max_sum > max)
		{
			max = post_max_sum;
			num_elems_max = j + 1;
		}
		//std::cout << "NE: " << j << " post sum: " << post_max_sum << "\n";
	}

	//std::cout << "MAX: " << max << " num elems: " << num_elems_max << "\n";
	// Put the selected into the partition with the root seq
	for(size_t j = 0; j < m_overlaps.size(); ++j)
	{
		if(j < num_elems_max)
			m_overlaps[j].partitionID = 0;
		else
			m_overlaps[j].partitionID = 1;
	}		

	double likelihood = calculateGroupedLikelihood();

	//std::cout << "MAX: " << max << " num elems: " << num_elems_max << "\n";
	// Put the selected into the partition with the root seq
	if(likelihood < initial_likelihood)
	{
		for(size_t j = 0; j < m_overlaps.size(); ++j)
			m_overlaps[j].partitionID = 0;
	}	
}

//
void MultiOverlap::partitionLI(double p_error)
{
	double initial_likelihood = calculateGroupedLikelihood();

	// Compute the likelihood of the alignment between the root sequence
	// and every member of the multioverlap
	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		double likelihood = 0.0f;
		double total_bases = 0;
		for(size_t j = 0; j < m_rootSeq.size(); ++j)
		{
			Pileup pileup = getSingletonPileup(j, i);
			if(pileup.getDepth() > 1)
			{
				DNADouble ap = pileup.calculateLikelihoodNoQuality(p_error);
				likelihood += ap.marginalize(0.25f);
				++total_bases;
			}
		}

		m_overlaps[i].score = likelihood / total_bases;
		//m_overlaps[i].partitionID = 1;
	}

	// Sort the overlaps by score
	std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::compareScore);

	std::cout << "Print score list\n";
	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		MOData& curr = m_overlaps[i];
		std::cout << i << "\t" << curr.ovr.id[1] << "\t" << curr.score << "\t" << curr.partitionID << "\n";
		//Initially set all elements to the "other" partition
		m_overlaps[i].partitionID = 1;
		
	}

	double max = -std::numeric_limits<double>::max();
	size_t num_elems_max = 0;

	double prev_likelihood = calculateGroupedLikelihood();

	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		m_overlaps[i].partitionID = 0;
		double likelihood = calculateGroupedLikelihood();
		//double improvement = likelihood - prev_likelihood;
		//if(improvement < 5)
		//	m_overlaps[i].partitionID = 0;
		//std::cout << "Elem: " << i << " prev: " << prev_likelihood << " likelihood: " << likelihood << " improvement: " << improvement << "\n";
		if(likelihood > max)
		{
			max = likelihood;
			num_elems_max = i + 1;
		}
		prev_likelihood = likelihood;
	}

	std::cout << "MAX: " << max << " NUM ELEMS: " << num_elems_max << "\n";
	// Put the selected into the partition with the root seq
	for(size_t j = 0; j < m_overlaps.size(); ++j)
	{
		if(max > initial_likelihood)
		{
			if(j < num_elems_max)
				m_overlaps[j].partitionID = 0;
			else
				m_overlaps[j].partitionID = 1;
		}
		else
		{
			m_overlaps[j].partitionID = 0;
		}
	}
	std::cout << "INIT: " << initial_likelihood << " MAX: " << max << " FINAL GL: " << calculateGroupedLikelihood() << "\n";	
}

//
bool MultiOverlap::partitionSL(double p_error, std::string dbg)
{
	(void)p_error;
	bool ret = false;
	std::string before;
	for(size_t j = 0; j < m_overlaps.size(); ++j)
	{
		if(m_overlaps[j].partitionID == 0)
			before.push_back('1');
		else
			before.push_back('0');
	}
	
	for(size_t j = 0; j < m_overlaps.size(); ++j)
	{
		m_overlaps[j].partitionID = 0;
	}


	std::cout << "\n\n **PartitionSL** \n\n";
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup p = getPileup(i);
		AlphaCount ac = p.getAlphaCount();

		// Estimate whether this pileup is a mixture between 
		// different source sequences
		// TODO: this is grossly oversimplified
		char sorted[ALPHABET_SIZE];
		ac.getSorted(sorted, ALPHABET_SIZE);
		//int most = ac.get(sorted[0]);
		int second = ac.get(sorted[1]);

		int cutoff = 3;
		if(second > cutoff)
		{
			// Generate mask
			char refBase = m_rootSeq[i];
			if(m_rootSeq[i] != dbg[i])
				ret = true;
			std::string mask;
			std::string str;
			for(size_t j = 0; j < m_overlaps.size(); ++j)
			{
				char b = getMODBase(m_overlaps[j], i);
				if(b == '\0' || b == refBase || (int)ac.get(refBase) < cutoff)
				{
					mask.push_back('1');
				}
				else
				{
					m_overlaps[j].partitionID = 1;
					mask.push_back('0');
				}

				if(b == '\0')
					str.push_back('-');
				else
					str.push_back(b);
			}
			std::cout << i << " " << refBase << dbg[i] << " MASK: " << mask << "\n";
			std::cout << i << " " << refBase << dbg[i] << " STR : " << str << "\n";
		}
	}

	std::string after;
	for(size_t j = 0; j < m_overlaps.size(); ++j)
	{
		if(m_overlaps[j].partitionID == 0)
			after.push_back('1');
		else
			after.push_back('0');
	}

	std::cout << "BEFORE: " << before << "\n";
	std::cout << " AFTER: " << after << "\n";
	return true;
}

//
bool MultiOverlap::partitionConflict(double p_error, std::string dbg)
{
	(void)p_error;
	std::cout << "\n\n **PartitionSL2** \n\n";
	std::vector<AlphaCount> acVec;
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup p = getPileup(i);
		AlphaCount ac = p.getAlphaCount();
		acVec.push_back(ac);
	
		char sorted[ALPHABET_SIZE];
		acVec[i].getSorted(sorted, ALPHABET_SIZE);
	}

	printf("ROOT\t*%s\t%d\n", m_rootSeq.c_str(), 0);
	printf("BASE\t*%s\t%d\n", dbg.c_str(), 0);

	for(size_t j = 0; j < m_overlaps.size(); ++j)
	{
		std::string cstr;
		std::string sc;
		std::string nc;

		int score = 0;
		int sum = 0;
		int wrong = 0;
		for(size_t i = 0; i < acVec.size(); ++i)
		{
			char rootBase = m_rootSeq[i];
			char sorted[ALPHABET_SIZE];
			acVec[i].getSorted(sorted, ALPHABET_SIZE);
			int second = acVec[i].get(sorted[1]);
			bool isConflict = second > 3;
			char b = getMODBase(m_overlaps[j], i);

			if(isConflict)
			{
				if(b != '\0' && (rootBase == sorted[0] || rootBase == sorted[1]))
				{
					if(b == rootBase)
						++score;
					else
						++wrong;
					++sum;
				}

				if(b == sorted[0] || b == sorted[1])
				{
					nc.push_back(b);
				}
				else
				{
					nc.push_back('N');
				}
			}
			
			if(b != '\0')
				cstr.push_back(b);	
			else
				cstr.push_back('N');
			
		}

		double frac;
		if(sum == 0)
			frac = 1.0f;
		else
			frac = (double)score/(double)sum;

		m_overlaps[j].score = frac;
		if((sum == 1 && wrong > 0) || (sum == 2 && wrong > 1) || wrong >= 2)
			m_overlaps[j].partitionID = 1;
		else
			m_overlaps[j].partitionID = 0;

		printf("CFT\t*%s\t%d,%lf\t%s\n", cstr.c_str(), (int)j, frac, nc.c_str());
	}

	/*
	std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::compareScore);
	for(size_t i = 0; i < m_overlaps.size() && i < 20; ++i)
	{
			m_overlaps[i].partitionID = 0;
	}
	*/
	return false;
}

//
std::string MultiOverlap::consensusConflict(double p_error)
{
	(void)p_error;

	const int CONFLICT_CUTOFF = 3;

	std::vector<AlphaCount> acVec;
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup p = getPileup(i);
		AlphaCount ac = p.getAlphaCount();
		acVec.push_back(ac);
	}

	for(size_t j = 0; j < m_overlaps.size(); ++j)
	{
		std::string cstr;
		std::string sc;
		std::string nc;

		int score = 0;
		int sum = 0;
		int wrong = 0;
		for(size_t i = 0; i < acVec.size(); ++i)
		{
			char rootBase = m_rootSeq[i];
			char sorted[ALPHABET_SIZE];
			acVec[i].getSorted(sorted, ALPHABET_SIZE);
			int second = acVec[i].get(sorted[1]);
			bool isConflict = second > CONFLICT_CUTOFF;
			char b = getMODBase(m_overlaps[j], i);

			if(isConflict)
			{
				if(b != '\0' && (rootBase == sorted[0] || rootBase == sorted[1]))
				{
					if(b == rootBase)
						++score;
					else
						++wrong;
					++sum;
				}

				if(b == sorted[0] || b == sorted[1])
				{
					nc.push_back(b);
				}
				else
				{
					nc.push_back('N');
				}
			}
			
			if(b != '\0')
				cstr.push_back(b);	
			else
				cstr.push_back('N');
		}

		double frac;
		if(sum == 0)
			frac = 1.0f;
		else
			frac = (double)score/(double)sum;

		m_overlaps[j].score = frac;
		if(sum > 0 && wrong > 0)
			m_overlaps[j].partitionID = 1;
		else
			m_overlaps[j].partitionID = 0;

		printf("CFT\t*%s\t%d,%lf\t%s\n", cstr.c_str(), (int)j, frac, nc.c_str());
	}

	std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::compareScore);
	for(size_t i = 0; i < m_overlaps.size() && i < 15; ++i)
	{
		m_overlaps[i].partitionID = 0;
	}

	std::string consensus = calculateConsensusFromPartition(p_error);
	for(size_t i = 0; i < acVec.size(); ++i)
	{
		char sorted[ALPHABET_SIZE];
		acVec[i].getSorted(sorted, ALPHABET_SIZE);
		int second = acVec[i].get(sorted[1]);
		bool isConflict = second > CONFLICT_CUTOFF;
		if(!isConflict)
		{
			consensus[i] = sorted[0];
		}
	}
	
	return consensus;
}

std::string MultiOverlap::consensusTemplate(const StringVec& templateVec)
{
	(void)templateVec;
	int maxScore = 0;
	int maxIdx = 0;
	for(size_t i = 0; i < templateVec.size(); ++i)
	{
		const std::string& tmpStr = templateVec[i];
		std::string m_str1;
		std::string m_str2;

		int score = 0;
		for(size_t j = 0; j < m_rootSeq.size(); ++j)
		{
			if(tmpStr[j] == m_rootSeq[j])
			{
				++score;
				m_str1.push_back('-');
				m_str2.push_back('-');
			}
			else
			{
				m_str1.push_back(m_rootSeq[j]);
				m_str2.push_back(tmpStr[j]);
			}
		}

		if(score > maxScore)
		{
			maxScore = score;
			maxIdx = i;
		}

		std::cout << "root\t" << m_str1 << "\n";
		std::cout << "tmpl\t" << m_str2 << "\t" << score << "\n";
	}
	std::cout << "selected " << maxIdx << " " << maxScore << "\n";
	return templateVec[maxIdx];
}



// Partition the MO using the best N elements only.
void MultiOverlap::partitionBest(double p_error, size_t n)
{
	// Compute the likelihood of the alignment between the root sequence
	// and every member of the multioverlap
	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		double likelihood = 0.0f;
		double total_bases = 0;
		double mismatches = 0;
		double overlap_len = 0;
		for(size_t j = 0; j < m_rootSeq.size(); ++j)
		{
			Pileup pileup = getSingletonPileup(j, i);
			if(pileup.getDepth() == 2)
			{
				++overlap_len;
				if(pileup.getBase(0) != pileup.getBase(1))
				{
					++mismatches;
				}
			}
			DNADouble ap = pileup.calculateLikelihoodNoQuality(p_error);
			likelihood += ap.marginalize(0.25f);
			total_bases += pileup.getDepth();
		}

		m_overlaps[i].score = likelihood / total_bases;
		m_overlaps[i].partitionID = 1;
		//std::cout << "ERS\t" << m_overlaps[i].score << "\t" << likelihood << "\n";
	}

	std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::compareScore);

	for(size_t i = 0; i < m_overlaps.size() && i < n; ++i)
	{
			m_overlaps[i].partitionID = 0;
	}
}

// Calculate the log probability that the string specified by mod originated from the string tmpStr
// base is the same frame of reference as m_rootSeq. The template string is assumed to be
// correct.
double MultiOverlap::calcTemplateProb(const std::string& tmpStr, double p_error, const MOData& mod) const
{
	double lp = 0.0f;
	double le = log(p_error);
	double lc = log(1 - p_error);
	for(size_t i = 0; i < tmpStr.size(); ++i)
	{
		char b = getMODBase(mod, i);
		if(b != '\0')
		{
			if(b == tmpStr[i])
				lp += lc;
			else
				lp += le;
		}	
	}
	return lp;
}

size_t MultiOverlap::countPartition(int id) const
{
	size_t count = 0;
	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		if(m_overlaps[i].partitionID == id)
			++count;
	}
	return count;
}

std::string MultiOverlap::calculateConsensusFromPartition(double p_error)
{
	std::string out;
	out.reserve(m_rootSeq.size());

	// require the best base call to be above this above to correct it
	double epsilon = 0.01;
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup p0;
		Pileup p1;
		getPartitionedPileup(i, p0, p1);
		DNADouble ap = p0.calculateLikelihoodNoQuality(p_error);
		char best_c;
		char curr_c = m_rootSeq[i];
		double max;
		ap.getMax(best_c, max);

		// Require the called value to be substantially better than the
		// current base
		if(best_c == curr_c || (best_c != curr_c && max - ap.get(curr_c) < epsilon))
		{
			out.push_back(curr_c);
		}
		else
		{
			out.push_back(best_c);
		}
	}
	return out;
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
		likelihood += ap.marginalize(0.25f);
	}
	return likelihood;
}

// Calculate the likelihood of the multioverlap given the groups defined by IntVector
// There are only two possible groups, 0 and 1. The root_sequence (index 0 in the vector) is assumed to belong
// to group 1
double MultiOverlap::calculateGroupedLikelihood() const
{
	double likelihood = 0.0f;
	for(size_t i = 0; i < m_rootSeq.size(); ++i)
	{
		Pileup g0;
		Pileup g1;
		getPartitionedPileup(i, g0, g1);
		if(g0.getDepth() > 0)
			likelihood += g0.calculateLikelihoodNoQuality(0.01).marginalize(0.25f);
		if(g1.getDepth() > 0)
			likelihood += g1.calculateLikelihoodNoQuality(0.01).marginalize(0.25f);
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

// Return the base in mod that matches the base at
// idx in the root seq. If mod does not overlap 
// the root at this position, returns '\0'
char MultiOverlap::getMODBase(const MOData& mod, int idx) const
{
	int trans_idx = idx - mod.offset;
	if(trans_idx >= 0 && size_t(trans_idx) < mod.seq.size())
		return mod.seq[trans_idx];
	else
		return '\0';
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

// Get the "stack" of bases that aligns to
// a single position of the root seq, including
// the root base only including the first numElems items
Pileup MultiOverlap::getPileup(int idx, int numElems) const
{
	assert((size_t)numElems < m_overlaps.size());
	Pileup p;
	p.add(m_rootSeq[idx]);

	for(int i = 0; i < numElems; ++i)
	{
		const MOData& curr = m_overlaps[i];
		// translate idx into the frame of the current sequence
		int trans_idx = idx - curr.offset;
		if(trans_idx >= 0 && (size_t)trans_idx < curr.seq.size())
		{
			p.add(curr.seq[trans_idx]);
		}
	}
	return p;
}


// Get the pileup between the root sequence and one of the sequences in the MO
Pileup MultiOverlap::getSingletonPileup(int base_idx, int ovr_idx) const
{
	assert(ovr_idx < (int)m_overlaps.size());

	Pileup p;
	p.add(m_rootSeq[base_idx]);

	const MOData& curr = m_overlaps[ovr_idx];
	// translate idx into the frame of the current sequence
	int trans_idx = base_idx - curr.offset;
	if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
	{
		p.add(curr.seq[trans_idx]);
	}
	return p;
}

// Fill in the pileup groups g0 and g1
void MultiOverlap::getPartitionedPileup(int idx, Pileup& g0, Pileup& g1) const
{
	g0.add(m_rootSeq[idx]);

	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		const MOData& curr = m_overlaps[i];
		// translate idx into the frame of the current sequence
		int trans_idx = idx - curr.offset;
		if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
		{
			if(curr.partitionID == 0)
				g0.add(curr.seq[trans_idx]);
			else
				g1.add(curr.seq[trans_idx]);
		}
	}
}


// Fill in the pileup groups g0 and g1
PileupVector MultiOverlap::getPartitionedPileup(int idx, int num_parts) const
{
	assert(false && "untested");
	PileupVector pv(num_parts);
	pv.push_back(Pileup());
	
	// the root always goes in partition 0
	pv.back().add(m_rootSeq[idx]);

	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		const MOData& curr = m_overlaps[i];
		// translate idx into the frame of the current sequence
		int trans_idx = idx - curr.offset;
		if(trans_idx >= 0 && size_t(trans_idx) < curr.seq.size())
		{
			assert(curr.partitionID < num_parts);
			pv[curr.partitionID].add(curr.seq[trans_idx]);
		}
	}
	return pv;
}


// Print the MultiOverlap to stdout
void MultiOverlap::print(int default_padding, int max_overhang)
{
	std::sort(m_overlaps.begin(), m_overlaps.end(), MOData::sortOffset);
	std::cout << "\nDrawing overlaps for read " << m_rootID << "\n";
	int root_len = int(m_rootSeq.size());
	
	// Print the root row at the bottom
	printRow(default_padding, max_overhang, root_len, 0, root_len, 0, 0.0f, m_rootSeq, m_rootID);

	for(size_t i = 0; i < m_overlaps.size(); ++i)
	{
		const MOData& curr = m_overlaps[i];
		int overlap_len = curr.ovr.match.getMaxOverlapLength();

		printRow(default_padding, max_overhang, root_len, 
		         curr.offset, overlap_len, curr.partitionID, 
				 curr.score, curr.seq, curr.ovr.id[1].c_str());
	}	
}

// Print the MultiOverlap groups specified by the IntVec to stdout
void MultiOverlap::printGroups()
{
	for(int i = 0; i <= 1; ++i)
	{
		MultiOverlap mo_group(m_rootID, m_rootSeq);
		for(size_t j = 0; j < m_overlaps.size(); ++j)
		{
			const MOData& curr = m_overlaps[j];
			if(curr.partitionID == i)
				mo_group.add(curr);
		}
		std::cout << "MO GROUP " << i << "\n";
		mo_group.print();
	}
}


// Print a single row of a multi-overlap to stdout
void MultiOverlap::printRow(int default_padding, int max_overhang, 
                            int root_len, int offset, int overlap_len, int pid, 
							double score, const std::string& seq, const std::string& id)
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
	printf("%s\t%d\t%d\t%lf\tID:%s\n", outstr.c_str(), overlap_len, pid, score, id.c_str());
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

