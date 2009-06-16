#include "FragmentDistribution.h"
#include <iostream>
#include <fstream>

void FragmentDistribution::readFromFile(std::string filename)
{
	m_totalCount = 0;
	std::ifstream reader(filename.c_str());
	int index;
	int count;
	while(reader >> index >> count)
	{
		m_counts.insert(std::make_pair(index, count));
		m_totalCount += count;
	}
	reader.close();
}

IntDist FragmentDistribution::convertToIntDist(double trim)
{
	int min, max;
	findPRange(trim, min, max);

	IntDist pdf(min, max);
	for(int idx = min; idx < max; ++idx)
	{
		pdf.addWeight(idx, getCount(idx));
	}
	pdf.normalize();
	return pdf;
}

void FragmentDistribution::findPRange(double p, int& min, int& max)
{
	double alpha = (1.0 - p) / 2;
	double sum = 0;
	for(IntIntMap::const_iterator iter = m_counts.begin(); iter != m_counts.end(); ++iter)
	{
		sum += getFreq(iter->first);
		if(sum > alpha)
		{
			min = iter->first;
			break;
		}
	}

	sum = 0.0f;
	for(IntIntMap::const_reverse_iterator iter = m_counts.rbegin(); iter != m_counts.rend(); ++iter)
	{
		sum += getFreq(iter->first);
		if(sum > alpha)
		{
			max = iter->first;
			break;
		}
	}
}

int FragmentDistribution::getCount(int index) const
{
	IntIntMap::const_iterator iter = m_counts.find(index);
	if(iter == m_counts.end())
	{
		return 0;
	}
	else
	{
		return iter->second;
	}
}

double FragmentDistribution::getFreq(int index) const
{
	return (double)getCount(index) / m_totalCount;
}
