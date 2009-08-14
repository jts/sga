#ifndef BWT_H
#define BWT_H

#include "STCommon.h"
#include "Occurance.h"

//
// BWT
//
class BWT
{
	public:
	
		// Constructors
		BWT(std::string t);
		BWT(const StringVector& sv);

		// Exact match
		void backwardSearch(std::string w);
		void getOverlaps(std::string w, int minOverlap);

		// Make all the cyclic rotations of a string
		static void makeCycles(SuffixString s, SuffixStringVector* outTable);

		// Print info about the BWT, including size
		void printInfo() const;

	private:

		void sortConstruct(int numStrings, SuffixStringVector* cycled);

		Occurance m_occurance;
		AlphaCount m_predCount;
		GSuffixVector m_suffixArray;
		
		std::string m_saStr;
		
		BWStr m_bwStr;
		int m_offset;
};

#endif
