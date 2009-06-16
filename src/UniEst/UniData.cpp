#include "UniData.h"
#include "StatsCommon.h"
#include <sstream>
#include <iostream>
#include <gsl/gsl_randist.h>

//
// Constructor
//
UniData::UniData(Contig c, int k) : m_contig(c),
							m_numReads(0),
							m_readCoverage(0),
							m_kmerCoverage(m_contig.getLength() - k + 1),
							m_depthUF(UF_NOCALL),
							m_overhangUF(UF_NOCALL)
{
	m_numOverhangingKmers[0] = 0;
	m_numOverhangingKmers[1] = 0;
	m_edgeMultiplicity[0] = 0;
	m_edgeMultiplicity[1] = 0;
}

//
// add coverage to the contig based on the alignment
//
void UniData::addCoverage(int k, const KAlignment& ka)
{
	++m_numReads;
	m_readCoverage += (ka.align_length);
	size_t start = ka.contig_start_pos;
	size_t end = start + (ka.align_length - k + 1);
	assert(start < m_kmerCoverage.size() && end <= m_kmerCoverage.size());

	for(size_t i = start; i < end; ++i)
		m_kmerCoverage[i]++;
}

//
// add "overhang" coverage - the set of k-mers linked to this contig by pairs
//
void UniData::addOverhangCoverage(int idx, int k, const KAlignment& ka)
{
	size_t start = ka.contig_start_pos;
	size_t end = start + (ka.align_length - k + 1);

	// Generate a string for each k-mer linked to the contig, based on its position
	for(size_t i = start; i < end; ++i)
	{
		std::stringstream ustr;
		ustr << ka.contig_id << ":" << i;
		m_pairedCoverage[idx].insert(ustr.str());
		m_numOverhangingKmers[idx]++;
	}
}

//
// Increment the number of edges in a particular direction
//
void UniData::addEdgeMultiplicity(int idx)
{
	assert(idx <= 1);
	m_edgeMultiplicity[idx]++;
}

//
// Calculate the expected number of unique overhanging kmers 
//
double UniData::getExpectedOverhangCoverage(const IntDist& fragDist, double mean) const
{
	// Calculate the expected kmer position distribution
	IntDist expectedDist = calculateExpectedOverhangDistribution(fragDist, getContigLen());
	//std::cout << expectedDist << std::endl;
	int start = expectedDist.getStart();
	int end = expectedDist.getEnd();
	int n = (int)(mean * (double)(end - start)); //m_numOverhangingKmers[idx];
	
	// For each position, calculate the probability of at least 1 kmer being seen = 1 - Bin(0, p, n)
	// The sum of these variables (each is a bernoulli random var) is the expected number
	double sum = 0.0f;
	for(int i = start; i < end; ++i)
	{
		double seen_prob = 1 - gsl_ran_binomial_pdf(0, expectedDist.getP(i), n);
		//double lp = n * log(1 - expectedDist.getP(i)
		//double seen_prob = 1 - (expectedDist.getP(start);

		//std::cout << "pos: " << i << " Prob of seeing: " << seen_prob << std::endl;
		sum += seen_prob;
	}

	//std::cout << "Expected: " << sum << " Observed: " << getOverhangCoverage(idx);
	return sum;
}

//
// Calculate the probability of seeing a kmer at distance d from the contig
// This is derived from the fragment distribution
// It ranges from 0 (the position immediately adjacent to the contig end) to fragDist.max (the furthest possible position)
// As each position in the contig is assumed to be equally likely to generate a PE fragment, this is
// just the (normalized) sum of the fragment distributions starting from each position
//
IntDist UniData::calculateExpectedOverhangDistribution(const IntDist& fragDist, int contigLength) const
{
	int fdStart = fragDist.getStart();
	int fdEnd = fragDist.getEnd();

	IntDist outDist(0, fdEnd - 1);

	int start = contigLength > fdEnd ? contigLength - fdEnd : 0;
	for(int i = start; i < contigLength; ++i)
	{
		for(int j = fdStart; j < fdEnd; ++j)
		{
			int end_pos = (i + j) - contigLength;
			if(end_pos >= 0)
			{
				outDist.addWeight(end_pos, fragDist.getP(j));
			}
		}
	}
	outDist.normalize();
	return outDist;
}

//
// Estimate whether the contig contained is unique based on the depth information
//
UniqueFlag UniData::estimateByDepth(double mean, double lrThreshold)
{
	int high_cn = 10;
	DoubleVec pVec;
	
	// Calculate the log-likelihood of the depth conditional on the copy number
	// P ( depth | cn ) = poisson(num reads, cn * mean * len)
	for(int cn = 1; cn < high_cn; ++cn)
	{
		double p1 = log_poisson(getNumReads(), cn * mean * getContigLen());
		//std::cout << mean_param << " " << data.getContigLen() << " " << data.getNumReads() << "\n";
		pVec.push_back(p1);
	}

	// Select the max value from the log-probability vector
	size_t best_cn = pVec.size();
	double best_p = 0.0f;
	for(size_t i = 0; i < pVec.size(); ++i)
	{
		if(best_cn == pVec.size() || pVec[i] > best_p)
		{
			best_cn = i + 1;
			best_p = pVec[i];
		}
	}

	// Calculate the likelihood ratio between the best and second best
	sort(pVec.begin(), pVec.end());
	double best = pVec[pVec.size() - 1];
	double second = pVec[pVec.size() - 2];
	double ratio = best - second;

	if(best_cn == 1 && ratio > lrThreshold)
		m_depthUF = UF_UNIQUE;
	else
		m_depthUF = UF_REPEAT;
	return m_depthUF;
}

//
// Estimate whether the contig is unique based on the overhanging paired data
//
UniqueFlag UniData::estimateByOverhang(const IntDist& fragDist, double mean, double error, double cutoff)
{
	double expected = getExpectedOverhangCoverage(fragDist, mean);
	double obs0 = (double)getOverhangCoverage(0) - (error * (double)getContigLen());
	double obs1 = (double)getOverhangCoverage(1) - (error * (double)getContigLen());;
	double ratio0 = (double)obs0 / expected;
	double ratio1 = (double)obs1 / expected;

	if(ratio0 <= cutoff && ratio1 <= cutoff)
		m_overhangUF = UF_UNIQUE;
	else
		m_overhangUF = UF_REPEAT;
	return m_overhangUF;
}

//
//
//
size_t UniData::getOverhangCoverage(int idx) const
{
	return m_pairedCoverage[idx].size();
}

//
//
//
double UniData::getMeanKmerCoverage() const
{
	double sum = 0.f;
	for(size_t i = 0; i < m_kmerCoverage.size(); ++i)
		sum += m_kmerCoverage[i];
	return sum / m_kmerCoverage.size();
}

//
//
//
void UniData::printStats() const
{
	std::cout << m_contig.getID() << "\t" 
				<< getContigLen() << "\t" 
				<< getNumReads() << "\t" 
				<< getReadDensity() << "\t" 
				<< m_edgeMultiplicity[0] << "\t"
				<< m_edgeMultiplicity[1] << "\t"
				<< m_depthUF << "\t"
				<< m_overhangUF << "\n";
}

//
//
//
void UniData::printKmerCoverage() const
{
	for(size_t i = 0; i < m_kmerCoverage.size(); ++i)
		std::cout << i << "\t" << m_kmerCoverage[i] << "\n";
}


