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
#ifndef SEARCHHISTORY_H
#define SEARCHHISTORY_H

#include "SearchTree.h"
#include "Util.h"

class SearchHistory
{
	public:
		
		//
		SearchHistory();

		//
		void add(int pos, char base);
		void normalize(bool doComplement);

		//
		static int countDifferences(const SearchHistory& a, const SearchHistory& b, int maxPos);

		//
		friend std::ostream& operator<<(std::ostream& out, const SearchHistory& hist);

	private:

		HistoryVector m_history;
};

#endif
