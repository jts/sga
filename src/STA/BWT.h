#ifndef BWT_H
#define BWT_H

#include "STCommon.h"


//
// Implementation of the occurance table B(a, i) which records the 
// number of times character a has been seen in B[0, i]
//
class OccuranceVec
{
	public:

		OccuranceVec() {}
		OccuranceVec(int s) {
			resize(s);
		}

		void resize(int s) {
			m_occ.resize(s);
		}

		int get(char a, int i) {
			return m_occ[i][a];
		}

		void set(char a, int i, int s) {
			m_occ[i][a] = s;
		}

		void inc(char a, int i) {
			m_occ[i][a]++;
		}

		void print()
		{
			for(size_t i = 0; i < m_occ.size(); i++)
			{
				for(std::map<char, int>::iterator iter = m_occ[i].begin(); iter != m_occ[i].end(); ++iter)
				{
					std::cout << "O(" << iter->first << "," << i << ") = " << iter->second << "\n";
				}
			}
		}

	private:
		std::vector<CharIntMap> m_occ;
};

//
// BWT
//
class BWT
{
	public:
	
		// Constructors
		BWT(std::string t);
		BWT(SAStringVector v);

		// Exact match
		void backwardSearch(std::string w);
		void getOverlaps(std::string w, int minOverlap);

		
		// Count the number of times c appears in string s[0, i]
		int countOcc(std::string s, int i, char c);

	private:

		void construct(int numStrings, SAStringVector& cycles);

		OccuranceVec m_occurances;
		CharIntMap m_predMap;
		IdentVector m_suffixArray;
		std::string m_saStr;
		std::string m_bwStr;
		std::string m_alphabet;
		int m_offset;
};

#endif
