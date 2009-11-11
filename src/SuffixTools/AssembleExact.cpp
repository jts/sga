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

		if(loop_count++ % 100 == 0)
			printf("[right] num active: %d loops: %d size: %zu\n", active_count, loop_count, pSeedList->size());
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
				leftExtend(minOverlap, curr, bv, pBWT, pRevBWT, pSAI);
			
			if(curr.isActive)
			{
				active = true;
				++active_count;
			}
		}

		if(loop_count++ % 100 == 0)
			printf("[left] num active: %d loops: %d size: %zu\n", active_count, loop_count, pSeedList->size());
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
	std::string& w = seed.seq;
	if(w.empty())
	{
		seed.isActive = false;
		return;
	}

	// The algorithm is as follows:
	// Starting from position 0 of w, we travel left to right finding the intervals in the BWT
	// and revBWT for each prefix i of w[0, i]. For each signficant prefix (length w[0,1] > minOverlap)
	// we determine the proper suffixes that match w[0,1]. We collect the spectrum of extensions for the 
	// all the suffixes that match a prefix of w (crucially this is constant time). If the extensions are unique
	// we prepend that to the seq, otherwise there is a branch at this sequence and we set isActive to false
	AlphaCount ext_counts;
	BWTIntervalPair ranges;
	BWTAlgorithms::initIntervalPair(ranges, w[0], pBWT, pRevBWT);

	for(size_t i = 1; i < w.length(); ++i)
	{
		/*
		std::string pre = w.substr(0, i);
		std::cout << "Ranges: " << ranges << "\n";
		std::cout << "pre: " << pre << "\n";
		*/
		BWTAlgorithms::updateBothR(ranges, w[i], pRevBWT);

		// Break if the string is no longer found
		if(!(ranges.interval[0].isValid() && ranges.interval[1].isValid())) 
			break;

		// ranges now holds the range in bwt/revbwt of the prefix w[0, i]
		if(i >= minOverlap)
		{
			// Map to the set of proper suffixes by calculating the range in 
			// pBWT that corresponds to w[0,1] + "$"
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothR(probe, '$', pRevBWT);
			
			// Calculate the spectrum of extensions within the range of valid suffixes
			// The difference between the occurance counts for the upper position and lower position
			// are all the possible extensions. The semantics of getOccDiff are sucht that the interval is closed
			// so that if you call getOccDiff(1,1) the result will be zeros. We want to include the 
			// endpoint values so we use upper + 1
			if(probe.interval[0].isValid())
			{
				assert(probe.interval[0].lower > 0);
				// Calculate the number of different extensions of the prefix S
				AlphaCount ac = pBWT->getOccDiff(probe.interval[0].lower - 1, probe.interval[0].upper);

				// Check if the $ count is set for this prefix
				size_t num_terminal = ac.get('$');
				
				if(num_terminal > 0)
				{
					// There may be multiple strings that terminate with this sequence
					// All these strings must have the same sequence and form a containment set
					// Since only the first string is (possibly) not contained (the others were marked in markInitIdentical)
					// we can just mark the first string in the terminal interval as being contained
					// Transform the position of the full-length suffix into it's read index
					size_t ri = bwtPos2ReadIdx(ranges.interval[0].lower, pBWT, pSAI);
					if(ri != seed.read_idx)
						bv[ri] = true;
				}
				ext_counts += ac;
				//printf("Proper suffix interval: %d, %d\n", (int)probe.interval[0].lower, (int)probe.interval[0].upper);
			}
		}
	}

	if(ext_counts.hasUniqueDNAChar())
	{
		char b = ext_counts.getUniqueDNAChar();
		seed.seq.insert(0, 1, b);
	}
	else
	{
		// Either there was a split in the possible extensions or there was no possible extension
		seed.isActive = false;
	}
}

// Extend the sequence by finding a consensus left-ward extension base
// and prepending it to seed.seq
void AssembleExact::rightExtend(const unsigned int minOverlap, ReadSeed& seed, bool_vec& bv, 
                                const BWT* pBWT, const BWT* pRevBWT, const SuffixArray* pSAI)
{
	//std::cout << "ID: " << seed.read_idx << " string: " << seed.seq << "\n";
	std::string& w = seed.seq;
	size_t l = w.length();
	int start = l - 1;

	if(w.empty())
	{
		seed.isActive = false;
		return;
	}

	AlphaCount ext_counts;
	BWTIntervalPair ranges;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);

	// The algorithm is as follows:
	// Starting from position l - 1 of w, we travel right to left finding the intervals in the BWT
	// and revBWT for each suffix i of w[i, l]. For each signficant suffix (length(w[i,l]) > minOverlap)
	// we determine the proper prefixes in the BWT that match the suffix. We collect the spectrum of extensions for the 
	// all the prefixes (crucially this is constant time). If the extensions are unique
	// we append that to the seq, otherwise there is a branch at this sequence and we set isActive to false
	for(int i = start - 1; i >= 0; --i)
	{
		// Update the ranges to give the range for the suffix w[i, l]
		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);

		std::string suf = w.substr(i);
		//std::cout << "Suf: " << suf << " range: " << ranges << "\n";

		// Break if the string is no longer found
		if(!(ranges.interval[0].isValid() && ranges.interval[1].isValid())) 
			break;

		// ranges now holds the range in bwt/revbwt of the prefix w[0, i]
		int o_len = l - i;
		if(o_len >= (int)minOverlap)
		{
			// Calculate which of the prefixes that match w[i, l] are terminal (they cannot be extended further)
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			
			// Calculate the spectrum of extensions within the range of valid suffixes
			// The difference between the occurance counts for the upper position and lower position
			// are all the possible extensions. The semantics of getOccDiff are sucht that the interval is closed
			// so that if you call getOccDiff(1,1) the result will be zeros. We want to include the 
			// endpoint values so we use upper + 1
			if(probe.interval[1].isValid())
			{
				assert(probe.interval[1].lower > 0);
				// Calculate the number of different extensions of the prefix S
				AlphaCount ac = pRevBWT->getOccDiff(probe.interval[1].lower - 1, probe.interval[1].upper);

				// Check if the $ count is set for this prefix
				size_t num_terminal = ac.get('$');
				
				if(num_terminal > 0)
				{
					// There may be multiple strings that terminate with this sequence
					// All these strings must have the same sequence and form a containment set
					// Since only the first string is (possibly) not contained (the others were marked in markInitIdentical)
					// we can just mark the first string in the terminal interval as being contained
					// Transform the position of the full-length suffix into it's read index
					//std::cout << "Suffix: " << w.substr(i) << " interval: " << ranges << "\n";
					size_t ri = bwtPos2ReadIdx(ranges.interval[0].lower, pBWT, pSAI);
					if(ri != seed.read_idx)
					{
						bv[ri] = true;
						//printf("marking string %s as terminal: %d\n", w.substr(0, i+1).c_str(), (int)ri);
					}
					
					static bool validation_flag = false;
					if(!validation_flag)
					{
						//printf("WARNING VALIDATION IS ENABLED IN LEFT EXTEND\n");
						validation_flag = true;
					}

					for(int64_t j = ranges.interval[0].lower; j <= ranges.interval[0].upper; ++j)
					{
						size_t ri = bwtPos2ReadIdx(j, pBWT, pSAI);
						if(ri != seed.read_idx)
							assert(bv[ri]);
					}	
				}
				ext_counts += ac;
				//printf("Proper suffix interval: %d, %d\n", (int)probe.interval[0].lower, (int)probe.interval[0].upper);
			}
		}
	}

	if(ext_counts.hasUniqueDNAChar())
	{
		char b = ext_counts.getUniqueDNAChar();
		seed.seq.append(1, b);
	}
	else
	{
		// Either there was a split in the possible extensions or there was no possible extension
		seed.isActive = false;
	}
}


#if 0
// Extend the sequence by finding a consensus right-ward extension base
// and appending it to seed.seq
void AssembleExact::rightExtend(const unsigned int minOverlap, ReadSeed& seed, const BWT* pBWT, const BWT* pRevBWT)
{
	std::string& w = seed.seq;
	size_t l = w.length();
	int start = w.length() - 1;

	//std::cout << "SEQ: " << w << "\n";
	if(w.empty())
	{
		seed.isActive[0] = false;
		seed.isActive[1] = false;
	}

	// This algorithm is fundamentally the same as leftExtend, see that function for a description
	AlphaCount ext_counts;
	BWTIntervalPair ranges;
	BWTAlgorithms::initIntervalPair(ranges, w[start], pBWT, pRevBWT);

	for(int i = start - 1; i >= 0; --i)
	{
		std::string suf = w.substr(i+1);

		BWTAlgorithms::updateBothL(ranges, w[i], pBWT);

		
		// Return if the string is no longer found
		if(!(ranges.interval[0].isValid() && ranges.interval[1].isValid())) 
			break;

		// ranges now holds the range in bwt/revbwt of the prefix w[0, i]
		int o_len = l - i;
		if(o_len >= (int)minOverlap)
		{
			// Map to the set of proper suffixes by calculating the range in 
			// pBWT that corresponds to w[0,1] + "$"
			BWTIntervalPair probe = ranges;
			BWTAlgorithms::updateBothL(probe, '$', pBWT);
			
			// Calculate the spectrum of extensions within the range of valid suffixes
			// The difference between the occurance counts for the upper position and lower position
			// are all the possible extensions. The semantics of getOccDiff are sucht that the interval is closed
			// so that if you call getOccDiff(1,1) the result will be zeros. We want to include the 
			// endpoint values so we use upper + 1
			if(probe.interval[1].isValid())
			{
				assert(probe.interval[1].lower > 0);
				AlphaCount ac = pRevBWT->getOccDiff(probe.interval[1].lower - 1, probe.interval[1].upper);
				ext_counts += ac;
				//printf("Proper suffix interval: %d, %d\n", (int)probe.interval[0].lower, (int)probe.interval[0].upper);
			}
		}
	}

	if(ext_counts.hasUniqueDNAChar())
	{
		char b = ext_counts.getUniqueDNAChar();
		seed.seq.append(1, b);
		//std::cout << "Appending " << b << "\n";
	}
	else
	{
		// Either there was a split in the possible extensions or there was no possible extension
		seed.isActive[1] = false;
	}	
}
#endif
