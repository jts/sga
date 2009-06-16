#ifndef UNIDATA_H
#define UNIDATA_H

#include "Util.h"
#include "Contig.h"
#include "IntDist.h"
#include <set>

//
// Class to hold data used to calculate whether a given contig is unique or not
//

//
// Typedef
//
typedef std::set<std::string> StringSet;


class UniData
{
	public:

		// Constructor
		UniData(Contig c, int k);

		// Add data
		void addCoverage(int k, const KAlignment& ka);
		void addOverhangCoverage(int idx, int k, const KAlignment& ka2);
		void addEdgeMultiplicity(int idx);

		// Get data
		Contig& getContig() { return m_contig; }
		const Contig& getContig() const { return m_contig; }

		size_t getContigLen() const { return m_contig.getLength(); }
		size_t getNumReads() const { return m_numReads; }
		size_t getReadCoverage() const { return m_readCoverage; }
		double getReadDensity() const { return (double)m_numReads / (double)m_contig.getLength(); } 
		size_t getOverhangCoverage(int idx) const;
		double getMeanKmerCoverage() const;
		size_t getEdgeMultiplicity(int idx) const { return m_edgeMultiplicity[idx]; }

		// Uniqueness estimate functions
		UniqueFlag estimateByDepth(double mean, double lrThreshold);
		UniqueFlag estimateByOverhang(const IntDist& fragDist, double mean, double error, double cutoff);

		// Calculate the number of expected overhanging kmers
		double getExpectedOverhangCoverage(const IntDist& fragDist, double mean) const;

		// Print
		void printStats() const;
		void printKmerCoverage() const;

	private:
		
		// Calculate the probability of seeing a kmer at distance d from the contig
		IntDist calculateExpectedOverhangDistribution(const IntDist& fragDist, int contigLength) const;

		// Member vars
		Contig m_contig;
		size_t m_numReads;
		size_t m_readCoverage;
		IntVec m_kmerCoverage;
		StringSet m_pairedCoverage[2];
		size_t m_numOverhangingKmers[2];
		size_t m_edgeMultiplicity[2];

		// These flags are used for information only, the actual unique is recorded in the contig
		// and determiend by a higher-level function
		UniqueFlag m_depthUF;
		UniqueFlag m_overhangUF;
};

#endif
