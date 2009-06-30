#ifndef PAIRSTREAMER_H
#define PAIRSTREAMER_H

#include "Util.h"
#include <iostream>
#include <fstream>

//
// Typedefs
//
typedef std::vector<AlignPair> AlignPairVec;

//
// Stream a pairs file from disc, one block of contigs at a time
//
class PairStreamer
{
	public:
		PairStreamer(std::string filename);
		~PairStreamer();
		AlignPairVec getBlock();

	private:
		AlignPair m_nextRecord;
		std::ifstream fileReader;
};

#endif
