//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// AssembleExact - Assembly algorithm for exact sequences using a BWT
//
#include "AssembleExact.h"

// Checks if the value in bv for a particular read seed is true
struct ContainChecker
{
	public:
		ContainChecker(const bool_vec* pBV) : m_pBV(pBV) {}
		bool operator()(const ReadSeed& rs)
		{
			return (*m_pBV)[rs.read_idx];
		}

	private:
		const bool_vec* m_pBV;
};

void AssembleExact::assemble(const unsigned int minOverlap, const BWT* pBWT, const BWT* pRevBWT, const ReadTable* pRT, const SuffixArray* pSAI)
{
	// Make a seed for each read
	RSList* pSeedList = new RSList;

	for(size_t i = 0; i < pRT->getCount(); ++i)
		pSeedList->push_back(ReadSeed(i, pRT->getRead(i).seq.toString()));

	// Initialize boolean array where bv[i] holds whether read i is contained within some other read
	std::vector<bool> bv(pRT->getCount(), false);

	std::cout << "Before initial contain: " << pSeedList->size() << "\n";

	// Find the initial containments
	markInitialIdentical(bv, pSeedList, pBWT, pSAI);
	
	// Remove the elements that are contained
	ContainChecker cc(&bv);
	pSeedList->remove_if(cc);
	std::cout << "After initial removal: " << pSeedList->size() << "\n";

	bool active = true;
	int loop_count = 0;
	int sample_rate = 2;
	while(active)
	{
		active = false;
		int active_count = 0;

		// Remove contained
		pSeedList->remove_if(cc);

		for(RSList::iterator iter = pSeedList->begin(); iter != pSeedList->end(); ++iter)
		{
			ReadSeed& curr = *iter;
			if(curr.isActive)
				leftExtend(minOverlap, curr, bv, pBWT, pRevBWT, pSAI);
			
			if(curr.isActive)
			{
				active = true;
				++active_count;
			}
		}

		if(loop_count++ % sample_rate == 0)
		{
			printf("[left] num active: %d loops: %d size: %zu\n", active_count, loop_count, pSeedList->size());
			sample_rate <<= 1;
		}
	}
	
	for(RSList::iterator iter = pSeedList->begin(); iter != pSeedList->end(); ++iter)
		iter->isActive = true;

	active = true;
	while(active)
	{
		active = false;
		int active_count = 0;

		// Remove contained
		pSeedList->remove_if(cc);

		for(RSList::iterator iter = pSeedList->begin(); iter != pSeedList->end(); ++iter)
		{
			ReadSeed& curr = *iter;
			if(curr.isActive)
				rightExtend(minOverlap, curr, bv, pBWT, pRevBWT, pSAI);
			
			if(curr.isActive)
			{
				active = true;
				++active_count;
			}
		}

		if(loop_count++ % sample_rate == 0)
		{
			printf("[right] num active: %d loops: %d size: %zu\n", active_count, loop_count, pSeedList->size());
			sample_rate <<= 1;
		}
	}
	
	std::ofstream out("assembled.fa");
	for(RSList::iterator iter = pSeedList->begin(); iter != pSeedList->end(); ++iter)
	{
		out << ">" << iter->read_idx << " " << iter->seq.length() << "\n";
		out << iter->seq << "\n";
	}
	out.close();

	delete pSeedList;
}

// Mark the initial containments, which are identical reads
void AssembleExact::markInitialIdentical(bool_vec& bv, const RSList* pList, const BWT* pBWT, const SuffixArray* pSAI)
{
	for(RSList::const_iterator iter = pList->begin(); iter != pList->end(); ++iter)
	{
		// Find the interval corresponding to the full sequence
		BWTInterval range = BWTAlgorithms::findInterval(pBWT, iter->seq);

		// Determine if the current node is the head of the block of substrings that are identical
		size_t s = bwtPos2ReadIdx(range.lower, pBWT, pSAI);

		// If this read is not the head of the block, mark its entry in bv to be true
		bv[iter->read_idx] = iter->read_idx != s;
	}
}

// If pos is the position in pBWT of a full-length suffix (the entire read), return its
// index into the suffix array index
// Invariant: pos is a full-length suffix, which is denoted by B[pos] == '$'
size_t AssembleExact::bwtPos2ReadIdx(int64_t pos, const BWT* pBWT, const SuffixArray* pSAI)
{
	assert(pBWT->getChar(pos) == '$');
	
	BWTInterval interval(pos, pos);
	BWTAlgorithms::updateInterval(interval, '$', pBWT);
	return pSAI->get(interval.lower).getID();
}

// Extend the sequence by finding a consensus left-ward extension base
// and prepending it to seed.seq
void AssembleExact::leftExtend(const unsigned int minOverlap, ReadSeed& seed, bool_vec& bv, 
                               const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI)
{
	std::string w = reverse(seed.seq.substr(0, seed.root_len));
	IntervalPairList* pContained = new IntervalPairList;
	if(w.empty())
	{
		seed.isActive = false;
		return;
	}

	AlphaCount ext_counts = collectExtensions(minOverlap, w, pRevBWT, pBWT, pContained);
	
	// Mark the heads of the intervals for the first entry of the pair (the intervals of pBWT) as being contained
	markHeadContained(bv, 1, seed.read_idx, pContained, pBWT, pSAI);
	
	if(ext_counts.hasUniqueDNAChar())
	{
		char b = ext_counts.getUniqueDNAChar();

		// There is a unique extension from the seed sequence to B
		// Ensure that the reverse is true and that the sequence we are extending to, extends to this one
		std::string joined_seq = b + seed.seq;
		std::string back_search = joined_seq.substr(0, seed.root_len);
		assert(back_search.length() == seed.root_len);

		// Get the count of bases back to the ending sequence of seed
		// We do this in the opposite direction of extension
		AlphaCount back_count = collectExtensions(minOverlap, back_search, pBWT, pRevBWT, NULL);

		if(back_count.hasUniqueDNAChar())
		{
			char r = back_count.getUniqueDNAChar();
			// Assert back is the character we are expecting
			assert(r == seed.seq[seed.root_len - 1]);
			seed.seq = joined_seq;
		}
		else
		{
			seed.isActive = false;
		}		
	}
	else
	{
		// Either there was a split in the possible extensions or there was no possible extension
		seed.isActive = false;
	}
	delete pContained;
}

// Extend the sequence by finding a consensus left-ward extension base
// and prepending it to seed.seq
void AssembleExact::rightExtend(const unsigned int minOverlap, ReadSeed& seed, bool_vec& bv, 
                                const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI)
{
	IntervalPairList* pContained = new IntervalPairList;
	std::string w = seed.seq.substr(seed.seq.length() - seed.root_len);

	if(seed.seq.empty())
	{
		seed.isActive = false;
		return;
	}

	AlphaCount ext_counts = collectExtensions(minOverlap, w, pBWT, pRevBWT, pContained);
	
	// Mark the heads of the intervals for the first entry of the pair (the intervals of pBWT) as being contained
	markHeadContained(bv, 0, seed.read_idx, pContained, pBWT, pSAI);
	
	if(ext_counts.hasUniqueDNAChar())
	{
		char b = ext_counts.getUniqueDNAChar();

		// There is a unique extension from the seed sequence to B
		// Ensure that the reverse is true and that the sequence we are extending to, extends to this one
		std::string joined_seq = seed.seq + b;
		std::string back_search = joined_seq.substr(joined_seq.length() - seed.root_len);

		// Get the count of bases back to the ending sequence of seed
		// We do this in the opposite direction of extension
		AlphaCount back_count = collectExtensions(minOverlap, reverse(back_search), pRevBWT, pBWT, NULL);

		if(back_count.hasUniqueDNAChar())
		{
			char r = back_count.getUniqueDNAChar();
			// Assert back is the character we are expecting
			assert(r == seed.seq[seed.seq.length() - seed.root_len]);
			seed.seq = joined_seq;
		}
		else
		{
			seed.isActive = false;
		}
	}
	else
	{
		// Either there was a split in the possible extensions or there was no possible extension
		seed.isActive = false;
	}

	delete pContained;
}

// Mark the head of each interval contained in pList as being contained
// Precondition: All entries in the interval that are not the head are 
// already marked as being contained. 
// Precondition: For each entry in the list, all the sequences for a given interval have an 
// identical proper prefix. When we perform the transformation bwtPos2ReadIdx we get ID
// of the head of a range of identical, full-length strings. Since all the non-head 
// identical strings were marked as contained in an earlier step, only the head is possibly
// active and we do not have to scan the entire range.
void AssembleExact::markHeadContained(bool_vec& bv, size_t int_idx, size_t read_idx, const IntervalPairList* pList, const BWT* pBWT, const SuffixArray* pSAI)
{
	for(IntervalPairList::const_iterator iter = pList->begin(); iter != pList->end(); ++iter)
	{
		const BWTInterval& interval = iter->interval[int_idx];
		size_t ri = bwtPos2ReadIdx(interval.lower, pBWT, pSAI);

		// Don't mark the index of the read that is currently being processed
		if(ri != read_idx)
			bv[ri] = true;
		
		static bool validation_flag = false;
		if(!validation_flag)
		{
			printf("WARNING VALIDATION IS ENABLED IN LEFT EXTEND\n");
			validation_flag = true;
		}

		for(int64_t j = interval.lower; j <= interval.upper; ++j)
		{
			size_t ri = bwtPos2ReadIdx(j, pBWT, pSAI);
			if(ri != read_idx)
				assert(bv[ri]);
		}
	}
}

// Return the count of all the possible one base extensions of the string w.
// This returns the number of times the suffix w[i, l]A, w[i, l]C, etc 
// appears in the FM-index for all i s.t. length(w[i, l]) >= minOverlap.
AlphaCount AssembleExact::collectExtensions(const unsigned int minOverlap, const std::string& w, 
                                            const BWT* pBWT, const BWT* pRevBWT, IntervalPairList* pIList)
{
	// The algorithm is as follows:
	// We perform a backward search on the FM-index of w.
	// For each signficant suffix (length w[i,l] >= minOverlap)
	// we determine the proper prefixes that match w[i,l]. For each proper prefix matching, 
	// we compute the number of extensions of A,C,G,T for those prefix.

	//std::cout << "ID: " << seed.read_idx << " string: " << seed.seq << "\n";
	
	AlphaCount ext_counts;
	BWTIntervalPair ranges;
	size_t l = w.length();
	int start = l - 1;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);

	for(int i = start - 1; i >= 0; --i)
	{
		// Compute the range of the suffix w[i, l]
		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);

		//std::cout << "Suf: " << w.substr(i) << " range: " << ranges << "\n";

		// Break if the suffix is no longer found
		if(!(ranges.interval[0].isValid() && ranges.interval[1].isValid())) 
			break;

		if((l - i) >= minOverlap)
		{
			// Calculate which of the prefixes that match w[i, l] are terminal
			// These are the proper prefixes (they are the start of a read)
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			
			// The probe interval contains the range of proper prefixes
			if(probe.interval[1].isValid())
			{
				assert(probe.interval[1].lower > 0);
				
				// The count for each extension is the difference between rank(B, upper) and rank(B, lower - 1)
				AlphaCount ac = pRevBWT->getOccDiff(probe.interval[1].lower - 1, probe.interval[1].upper);

				// If any of the prefixes terminated with the $ character, it indicates 
				// that read it corresponds to is contained within the current sequence
				// The interval that contains these proper suffixes is returned in pIList
				// if pIList is not NULL
				size_t num_terminal = ac.get('$');
				
				if(num_terminal > 0 && pIList != NULL)
					pIList->push_back(ranges);
				
				ext_counts += ac;
				//printf("Proper suffix interval: %d, %d\n", (int)probe.interval[0].lower, (int)probe.interval[0].upper);
			}
		}
	}
	return ext_counts;
}
