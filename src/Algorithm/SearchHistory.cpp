//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SearchHistory - This class stores the sequence
// of discordant bases that were used while searching
// for a string in the FM-index
//
#include "SearchHistory.h"
#include "Util.h"
#include <iterator>
#include <algorithm>
//
SearchHistory::SearchHistory()
{

}

// Normalize the search history so that the positions are in ascending (left to right)
// order and the bases are wrt. the original strand of the forward sequence
// This is required to compare two SearchHistories that represent alignments to different 
// strands
void SearchHistory::normalize(bool doComplement)
{
	std::sort(m_history.begin(), m_history.end(), SearchHistoryItem::sortPos);

	if(doComplement)
	{
		for(size_t i = 0; i < m_history.size(); ++i)
		{
			m_history[i].base = complement(m_history[i].base);	
		}
	}
}

//
void SearchHistory::add(int pos, char base)
{
	SearchHistoryItem item(pos, base);
	m_history.push_back(item);
}

//
int SearchHistory::countDifferences(const SearchHistory& a, const SearchHistory& b, int maxPos)
{
	size_t na = a.m_history.size();
	size_t nb = b.m_history.size();

	size_t i = 0;
	size_t j = 0;
	int count = 0;

	while(i < na && j < nb)
	{
		const SearchHistoryItem& itemA = a.m_history[i];
		const SearchHistoryItem& itemB = b.m_history[j];
		
		if(itemA.pos > maxPos || itemB.pos > maxPos)
			break;

		if(itemA.pos == itemB.pos)
		{
			if(itemA.base != itemB.base)
				++count;
			++i;
			++j;
		}
		else if(itemA.pos < itemB.pos)
		{
			++count;
			++i;
		}
		else if(itemB.pos < itemA.pos)
		{
			++count;
			++j;
		}
		else
		{
			assert(false);
		}
	}
	
	// Count the remaining elements in A and B that are less than maxpos
	while(i < na && a.m_history[i].pos <= maxPos)
	{
		++count;
		++i;
	}

	// Count the remaining elements in A and B that are less than maxpos
	while(j < nb && b.m_history[j].pos <= maxPos)
	{
		++count;
		++j;
	}

	return count;
}

//
std::ostream& operator<<(std::ostream& out, const SearchHistory& hist)
{
	std::copy(hist.m_history.begin(), hist.m_history.end(), 
	          std::ostream_iterator<SearchHistoryItem>(out, "\t"));
	return out;
}

