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
#include <iterator>

SearchHistory::SearchHistory()
{

}

void SearchHistory::add(int pos, char base)
{
	HistoryItem item = {pos, base};
	m_history.push_back(item);
}

std::ostream& operator<<(std::ostream& out, const SearchHistory& hist)
{
	std::copy(hist.m_history.begin(), hist.m_history.end(), 
	          std::ostream_iterator<SearchHistory::HistoryItem>(out, "\t"));
	return out;
}

