//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SuffixCompare - Suffix Comparator class
// Designed to be passed to bucket/histogram sort and std::sort
//
#ifndef SUFFIXCOMPARE_H
#define SUFFIXCOMPARE_H
#include "STCommon.h"
#include "ReadTable.h"

class LCPArray;

// Simple struct to sort suffixes using their ids by looking
// up substrings in the table via pRT
struct SuffixCompare
{
	public:
		SuffixCompare(const ReadTable* pRT);
		SuffixCompare(const ReadTable* pRT, int offset, int bucket_len);
		SuffixCompare(const SuffixCompare& other);
		~SuffixCompare();

		// Comparator function
		bool operator()(SAElem x, SAElem y) const; 

		// Bucket function
		int operator()(SAElem x) const;
		
		// Calculate the number of possible suffixes
		int calcNumSuffixes(int maxLen) const;

		// Return the number of buckets needed for bucket sort
		int getNumBuckets() const;

		// Return the bucket offset
		size_t getBucketOffset() const { return m_bucketOffset; }

		// Set the offset
		void setBucketDepth(int depth);

		// Return true if a bucket is degenerate ie it cannot be subdivided any further
		bool isBucketDegenerate(int index) const;

		// Calculate the number of suffixes that precede the first instance of b for a 
		// given maximum suffix length
		int numPredSuffixes(char b, int maxLen) const;

	private:

		inline uint8_t getRank(char b) const
		{
			return m_rankLUT[static_cast<uint8_t>(b)];
		}
		SuffixCompare() {}
		
		const ReadTable* m_pRT;
		size_t m_bucketOffset;
		size_t m_bucketLen;

		static const uint8_t m_rankLUT[256];
};

#endif
