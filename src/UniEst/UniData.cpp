#include "UniData.h"
#include <sstream>
#include <iostream>

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

void UniData::addOverhangCoverage(int idx, int k, const KAlignment& ka)
{
	size_t start = ka.contig_start_pos;
	size_t end = start + (ka.align_length - k + 1);

	for(size_t i = start; i < end; ++i)
	{
		std::stringstream ustr;
		ustr << ka.contig_id << ":" << i;
		m_pairedCoverage[idx].insert(ustr.str());
		//std::cout << "posstr: " << ustr.str() << std::endl;
	}
}

void UniData::printKmerCoverage() const
{
	for(size_t i = 0; i < m_kmerCoverage.size(); ++i)
		std::cout << i << "\t" << m_kmerCoverage[i] << "\n";
}

double UniData::getMeanKmerCoverage() const
{
	double sum = 0.f;
	for(size_t i = 0; i < m_kmerCoverage.size(); ++i)
		sum += m_kmerCoverage[i];
	return sum / m_kmerCoverage.size();
}

