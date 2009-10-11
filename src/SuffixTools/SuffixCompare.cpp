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
	const char* suffix_x = rx.seq.getSuffix(x.getPos());
	const char* suffix_y = ry.seq.getSuffix(y.getPos());
	
	int cmp = strcmp(suffix_x, suffix_y);

	// If the suffixes are identical all the way to the last char, break ties by id
	if(cmp != 0)
		return cmp < 0;
	else
		return rx.id < ry.id;
}

// Get the bucket for a particular SAElem
int SuffixCompare::operator()(SAElem x) const
{
	//std::cout << "Finding bucket for " << x << "\n";
	const DNAString& read = m_pRT->getRead(x.getID()).seq;

	size_t suffix_start = x.getPos() + m_bucketOffset;
	const char* suffix = read.getSuffix(suffix_start);
	size_t suffix_len = read.getSuffixLength(suffix_start);

	size_t stop = std::min(m_bucketLen, suffix_len);

	int rank = 0;
	for(size_t i = 0; i < stop; ++i)
	{
		char b = suffix[i];
		rank += numPredSuffixes(b, m_bucketLen - i);
	}

	return rank;
}

//
void SuffixCompare::printElem(SAElem& x) const
{
	const DNAString& read = m_pRT->getRead(x.getID()).seq;
	std::cout << read.getSuffixString(x.getPos()) << "\n";
}

//
void SuffixCompare::setBucketDepth(int depth)
{
	m_bucketOffset = m_bucketLen * depth;
}

//
int SuffixCompare::getNextSortingIndex() const
{
	return m_bucketOffset + 1;
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

