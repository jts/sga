//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapQuality - Probabilitic inference of the quality of a set of overlaps
//
#include "OverlapQuality.h"
#include <math.h>
#include <limits>

OverlapQuality::OverlapQuality(const SeqItem& base_read, const OverlapVector& overlaps, const ReadTable* pRT)
{
	double p_err = 0.01f;

	m_baseSeq = base_read.seq.toString();
	m_baseQual.reserve(m_baseSeq.length());
	for(size_t i = 0; i < m_baseSeq.length(); ++i)
		m_baseQual[i] = p_err;

	m_cSeqs.reserve(overlaps.size());
	m_cQuals.reserve(overlaps.size());

	for(size_t i = 0; i < overlaps.size(); ++i)
	{
		Overlap curr = overlaps[i];

		// If the first element of the overlap is not the read
		// swap the elements
		if(curr.id[0] != base_read.id)
			curr.swap();

		assert(curr.id[0] == base_read.id);
	
		std::string otherSeq = pRT->getRead(curr.id[1]).seq.toString();

		// Make the other sequence in the same frame as the root
		if(curr.match.isRC())
		{
			otherSeq = reverseComplement(otherSeq);
			curr.match.canonize();
		}

		double_vec r;
		std::string c;
		c.reserve(m_baseSeq.length());
		r.reserve(m_baseSeq.length());

		for(size_t j = 0; j < m_baseSeq.length(); ++j)
		{
			int transformed = curr.match.translate(j);
			if(transformed < 0 || transformed >= (int)otherSeq.length())
			{
				c.append(1, '.');
				r.push_back(1.0f);
			}
			else
			{
				c.append(1, otherSeq[transformed]);
				r.push_back(p_err);
			}
		}
		m_cSeqs.push_back(c);
		m_cQuals.push_back(r);
	}
}

// Calculate the likelihood of the data given that the true sequence is x
double OverlapQuality::calcLikelihood(const::std::string& x) const
{
	double sum = 0.0;

	// Calculate the contribution of the base sequence first
	sum = calcAlignL(x, m_baseSeq, m_baseQual);

	// Calculate the contribution from the aligning reads
	for(size_t i = 0; i < m_cSeqs.size(); ++i)
	{
		double al = calcAlignL(x, m_cSeqs[i], m_cQuals[i]);
		sum += al;
	}
	return sum;
}

double OverlapQuality::cluster()
{
	printf("[cluster] base seq: %s\n", m_baseSeq.c_str());
	// Set up the initial cluster index values
	uint_vec idx_vec(m_cSeqs.size(), 0);

	// Set up the initial template
	StringVector t_vec;

	std::string initial_template = posteriorDecode(0, idx_vec);
	t_vec.push_back(initial_template);

	// Calculate the initial likelihood of all the data
	double init_l = calculateClusterLikelihood(t_vec, idx_vec);

	double min_s = std::numeric_limits<double>::max();
	size_t min_idx = 0;
	
	for(size_t i = 0; i < m_cSeqs.size(); ++i) 
	{
		size_t k = idx_vec[i];
		
		double al = calcAlignL(t_vec[k], m_cSeqs[i], m_cQuals[i]);
		double tl = calcTrueL(m_cSeqs[i], m_cQuals[i]);

		double score = al - tl;
		if(score < min_s)
		{
			min_s = score;
			min_idx = i;
		}
	}

	idx_vec[min_idx] = t_vec.size();
	std::string temp = posteriorDecode(t_vec.size(), idx_vec);

	printf("Choosing new template: score %lf idx %zu\n", min_s, min_idx);
	printf("Temp: %s\n", temp.c_str());
	t_vec.push_back(temp);

	for(size_t i = 0; i < t_vec.size(); ++i)
		printf("T(%zu): %s\n", i, t_vec[i].c_str());

	// Reassign the sequences to the templates
	for(size_t i = 0; i < m_cSeqs.size(); ++i)
	{
		double max_l = -std::numeric_limits<double>::max();
		size_t max_k = 0;
		for(size_t k = 0; k < t_vec.size(); ++k)
		{
			double al = calcAlignL(t_vec[k], m_cSeqs[i], m_cQuals[i]);
			if(al > max_l)
			{
				max_l = al;
				max_k = k;
			}
			printf("%zu al: %lf\n", k, al);
		}

		idx_vec[i] = max_k;
		printf("C(%u): %s\n", idx_vec[i], getMaskedSeq(m_cSeqs[i]).c_str());
	}

	// recalculate the posterior template
	t_vec[1] = posteriorDecode(1, idx_vec);

	// Calculate the likelihood given the new assignments
	double final_l = calculateClusterLikelihood(t_vec, idx_vec);
	double improve = final_l - init_l;
	printf("init l: %lf final l: %lf improve: %lf\n", init_l, final_l, improve);
	return improve;
}

std::string OverlapQuality::posteriorDecode(size_t k, uint_vec& idx_vec)
{
	std::string out;
	out.reserve(m_baseSeq.length());

	// Loop over each position of the base sequence and calculate the 
	// posterior probability of each possible base given all the sequences
	// that are in set k (marked by idx_vec)
	for(size_t j = 0; j < m_baseSeq.length(); ++j)
	{
		int count = 0;

		double posterior[DNA_ALPHABET_SIZE];
		for(size_t b = 0; b < DNA_ALPHABET_SIZE; ++b)
			posterior[b] = 0.0f;

		for(size_t i = 0; i < m_cSeqs.size(); ++i)
		{
			if(idx_vec[i] != k)
				continue;
		
			for(size_t b = 0; b < DNA_ALPHABET_SIZE; ++b)
			{
				if(m_cSeqs[i][j] == ALPHABET[b])
					posterior[b] += log(1.0 - m_cQuals[i][j]);
				else
					posterior[b] += log(m_cQuals[i][j]);
			}
		
			if(m_cSeqs[i][j] != '.')
				++count;

		}

		// Make sure there was some valid data at this point, otherwise use the base seq
		if(count == 0)
			out.append(1, m_baseSeq[j]);
		else
		{
			// Use the maximum posterior value, assuming each one is equally likely
			size_t max_b = 0;
			double max_p = -std::numeric_limits<double>::max();

			for(size_t b = 0; b < DNA_ALPHABET_SIZE; ++b)
			{
				if(posterior[b] > max_p)
				{
					max_p = posterior[b];
					max_b = b;
				}
			}
			//printf("maxp: %lf b: %c\n", max_p, ALPHABET[max_b]);
			out.append(1, ALPHABET[max_b]);
		}
	}
	return out;
}

// Calculate the total likelihood of the data given the clustering
double OverlapQuality::calculateClusterLikelihood(const StringVector& t_vec, const uint_vec& idx_vec)
{
	double sum = 0.0f;
	for(size_t i = 0; i < m_cSeqs.size(); ++i) 
	{
		size_t k = idx_vec[i];
		sum += calcAlignL(t_vec[k], m_cSeqs[i], m_cQuals[i]);
	}
	return sum;
}

double OverlapQuality::calcAlignL(const std::string& x, const std::string& c, const double_vec& r) const
{
	double sum = 0.0;
	for(size_t j = 0; j < x.length(); ++j)
		if(x[j] != '.')
			sum += x[j] == c[j] ? log(1.0f - r[j]) : log(r[j]);
	return sum;
}

double OverlapQuality::calcTrueL(const std::string& c, const double_vec& r) const
{
	double sum = 0.0;
	for(size_t j = 0; j < r.size(); ++j)
	{
		if(c[j] != '.')
			sum += log(1.0f - r[j]);
	}
	return sum;
}

// mask out the positions that match the base sequence
std::string OverlapQuality::getMaskedSeq(const std::string& str) const
{
	std::string out(str);
	for(size_t i = 0; i < str.length(); ++i)
	{
		if(str[i] != m_baseSeq[i] ||str[i] == '.')
			out[i] = str[i];
		else
			out[i] = ',';
	}
	return out;
}


