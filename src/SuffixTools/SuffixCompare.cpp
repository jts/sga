//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SuffixCompare 
//
#include "SuffixCompare.h"

const uint8_t SuffixCompare::m_rankLUT[256] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,0,
	0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

//
SuffixCompare::SuffixCompare(const ReadTable* pRT) : m_pRT(pRT)
{
	m_bucketOffset = 0;
	m_bucketLen = 6;
}

//
SuffixCompare::SuffixCompare(const ReadTable* pRT, int offset, int bucket_len) : m_pRT(pRT),
                                                                                 m_bucketOffset(offset),
                                                                                 m_bucketLen(bucket_len)
{
	if(m_bucketOffset % m_bucketLen != 0)
	{
		std::cerr << "Error, bucket offset must be a multiple of bucket length\n";
		assert(false);
	}
}

//
SuffixCompare::SuffixCompare(const SuffixCompare& other)
{
	m_pRT = other.m_pRT;
	m_bucketOffset = other.m_bucketOffset;
	m_bucketLen = other.m_bucketLen;
}
//
SuffixCompare::~SuffixCompare()
{
}

// Compare two suffixes
bool SuffixCompare::operator()(SAElem x, SAElem y) const
{ 
	const SeqItem& rx = m_pRT->getRead(x.getID());
	const SeqItem& ry = m_pRT->getRead(y.getID());
	const std::string& sx = rx.seq;
	const std::string& sy = ry.seq;
	
	std::string sfx = sx.substr(x.getPos()) + "$";
	std::string sfy = sy.substr(y.getPos()) + "$";
	int cmp = sfx.compare(sfy);
	if(cmp == 0)
		return rx.id < ry.id;
	else
		return cmp < 0;
}

// Get the bucket for a particular SAElem
int SuffixCompare::operator()(SAElem x) const
{
	//std::cout << "Finding bucket for " << x << "\n";
	const std::string& r = m_pRT->getRead(x.getID()).seq;

	size_t position = x.getPos() + m_bucketOffset;
	std::string sfx;
	if(position >= r.length())
		sfx = "$";
	else
		sfx = r.substr(position) + "$";

	size_t stop = std::min(m_bucketLen, sfx.length());
	int rank = 0;
	for(size_t i = 0; i < stop; ++i)
	{
		char b = sfx[i];
		rank += numPredSuffixes(b, m_bucketLen - i);
	}
	return rank;
}

//
void SuffixCompare::setBucketDepth(int depth)
{
	m_bucketOffset = m_bucketLen * depth;
}

//
int SuffixCompare::calcNumSuffixes(int maxLen) const
{
	int r = 0;
	for(int i = 0; i <= maxLen; ++i)
		r += (1 << 2*i);
	return r;
}

bool SuffixCompare::isBucketDegenerate(int index) const
{
	if(index == 0)
		return true;
	else
	{
		for(int i = m_bucketLen; i >= 2; --i)
		{
			--index;
			int block_size = calcNumSuffixes(i - 1);
			if(index % block_size == 0)
				return true;
			else
			{
				int block = index / block_size;
				index -= block * block_size;
			}
		}
	}
	return false;
}

//
int SuffixCompare::getNumBuckets() const
{
	return calcNumSuffixes(m_bucketLen);
}

// Returns the number of suffixes that are less than the base b for the given max length
int SuffixCompare::numPredSuffixes(char b, int maxLen) const
{
	// base case
	int rb = getRank(b);
	if(rb == 0)
		return 0;
	int block_size = calcNumSuffixes(maxLen - 1);
	return block_size * (rb - 1) + 1;
}

