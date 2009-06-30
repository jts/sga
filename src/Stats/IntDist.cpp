#include "IntDist.h"
#include <iostream>

double IntDist::ERROR_LIMIT = 1e-16;

//
//
//
IntDist::IntDist(int min, int max) : m_start(min), m_normalized(true)
{
	assert(max >= min);
	m_values.resize(max - min + 1, 0.0f);
}

//
// Set the probabilty by giving a weight first (normalize must be called later)
//
void IntDist::addWeight(int pos, double w)
{
	size_t idx = pos2Idx(pos);
	assert(idx < m_values.size());
	m_values[idx] += w;
	m_normalized = false;
}

//
// Calculate the expected value
// 
double IntDist::expectedValue() const
{
	double sum = 0.0f;
	for(int i = getStart(); i <= getEnd(); ++i)
		sum += (double)(i) * getP(i);
	return sum;
}
//
// Normalize the pdf by converting weights to probabilities
//
void IntDist::normalize()
{
	double sum = 0.0f;
	for(size_t idx = 0; idx < m_values.size(); ++idx)
		sum += m_values[idx];
	
	if(sum > ERROR_LIMIT)
	{
		for(size_t idx = 0; idx < m_values.size(); ++idx)
			m_values[idx] /= sum;
	}
	else
	{
		// The distribution is degenerate, set all values to zero
		for(size_t idx = 0; idx < m_values.size(); ++idx)
			m_values[idx] = 0;
	}	

	m_normalized = true;
}

//
// Set the probability directly
//
void IntDist::setP(int pos, double p)
{
	size_t idx = pos2Idx(pos);
	assert(idx < m_values.size());
	m_values[idx] = p;
}

//
// Get the probability at position pos
// The distribution must be normalized
//
double IntDist::getP(int pos) const
{
	assert(m_normalized);
	return getWeight(pos);
}

//
// Get the weight at position pos
// The distribution is not necessarily normalized
//
double IntDist::getWeight(int pos) const
{
	size_t idx = pos2Idx(pos);
	if(idx == m_values.size()) // out-of-range position
		return 0.0f;
	return m_values[idx];
}


//
// Convert the position to the index into the vector
//
size_t IntDist::pos2Idx(int pos) const
{
	size_t idx = pos - m_start;
	if(pos < m_start || idx >= m_values.size())
	{
		return m_values.size();
	}
	return idx;
}

//
// Write out the pdf
//
std::ostream& operator<<(std::ostream& out, const IntDist& pdf)
{
	for(size_t idx = 0; idx < pdf.m_values.size(); ++idx)
	{
		out << (int)idx + pdf.m_start << "\t" << pdf.m_values[idx] << "\n";
	}
	return out;
}

