#ifndef UNIDATA_H
#define UNIDATA_H

#include "Util.h"
#include "Contig.h"
#include <set>

//
// Typedef
//
typedef std::set<std::string> StringSet;


class UniData
{
	public:
		UniData(Contig c, int k) : m_contig(c),
									m_numReads(0),
									m_readCoverage(0),
									m_kmerCoverage(m_contig.getLength() - k + 1) {}

		//
		void addCoverage(int k, const KAlignment& ka);
		void addOverhangCoverage(int idx, int k, const KAlignment& ka2);
		
		//
		Contig& getContig() { return m_contig; }
		size_t getContigLen() const { return m_contig.getLength(); }
		size_t getNumReads() const { return m_numReads; }
		size_t getReadCoverage() const { return m_readCoverage; }
		double getReadDensity() const { return (double)m_numReads / (double)m_contig.getLength(); } 

		//
		double getMeanKmerCoverage() const;
		void printKmerCoverage() const;


	private:
		Contig m_contig;
		size_t m_numReads;
		size_t m_readCoverage;
		IntVec m_kmerCoverage;
		StringSet m_pairedCoverage[2];
};

#endif
