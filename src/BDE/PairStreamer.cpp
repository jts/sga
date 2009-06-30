#include "PairStreamer.h"

//
//
//
PairStreamer::PairStreamer(std::string filename)
{
	fileReader.open(filename.c_str());
	assert(fileReader.is_open());

	// Prime the read by taking in the first record
	fileReader >> m_nextRecord;
}

//
// Read in a block of records (records that all have the same contig_id for the
// first alignment). 
//
AlignPairVec PairStreamer::getBlock()
{
	AlignPairVec out;

	// If the previous block ended by hitting the end of the file, the EOF
	// flag will be set, return an empty array if so
	if(fileReader.eof())
	{
		return out;
	}

	// There are still alignments to read, push the first record
	out.push_back(m_nextRecord);

	AlignPair currRecord;
	while(fileReader >> currRecord)
	{
		// The read was successful 
		if(currRecord.aligns[0].contig_id == m_nextRecord.aligns[0].contig_id)
		{
			// Still in the block
			out.push_back(currRecord);
		}
		else
		{
			// block ended
			m_nextRecord = currRecord;
			break;
		}
	}
	return out;
}

//
//
//
PairStreamer::~PairStreamer()
{
	fileReader.close();
}


