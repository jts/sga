//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// BWTDiskConstruction - Implementation of the
// BWT disk construction algorithm of Ferragina,
// Gagie, Manzini. See Lightweight Data Indexing
// and Compression in External Memory
//
#include "BWTDiskConstruction.h"
#include "Util.h"

#define WORKING_FILENAME "partialbwt."

struct MergeItem
{
	size_t start_index;
	size_t end_index;
	std::string filename;

	friend std::ostream& operator<<(std::ostream& out, const MergeItem& item)
	{
		out << "[" << item.start_index << "," << item.end_index << "] " << item.filename;
		return out;
	}
};

typedef std::vector<MergeItem> MergeVector;

// Here we create M bwts for the input
// file where M depends on the amount of memory
// available. After these have been created
// we merge the bwts in pairs
BWT* buildBWTDisk(const std::string& filename, const std::string& bwt_extension)
{
	size_t MAX_READS_PER_GROUP = 77;

	SeqReader reader(filename);
	SeqRecord record;

	ReadTable* pCurrRT = new ReadTable;

	int groupID = 0;
	size_t numReadTotal = 0;

	MergeVector mergeVector;
	MergeItem mergeItem;
	mergeItem.start_index = 0;

	// Compute the initial BWTs
	
	bool done = false;
	while(!done)
	{
		done = !reader.get(record);

		if(!done)
		{
			// the read is valid
			pCurrRT->addRead(record.toSeqItem());
			++numReadTotal;
		}

		if(pCurrRT->getCount() >= MAX_READS_PER_GROUP || done)
		{
			if(pCurrRT->getCount() > 0)
			{
				// Compute the SA and BWT for this group
				SuffixArray* pSA = new SuffixArray(pCurrRT);
				BWT* pBWT = new BWT(pSA, pCurrRT);

				// Write the BWT to disk				
				std::stringstream ss;
				ss << WORKING_FILENAME << groupID++ << bwt_extension;
				std::string temp_filename = ss.str();
				pBWT->write(temp_filename);
				delete pBWT;
				delete pSA;
				delete pCurrRT;

				// Save the item for later merging
				mergeItem.end_index = numReadTotal - 1;
				mergeItem.filename = temp_filename;
				mergeVector.push_back(mergeItem);

				// Reset the data
				mergeItem.start_index = numReadTotal;
				if(!done)
					pCurrRT = new ReadTable;
			}
			else
			{
				delete pCurrRT;
			}
		}
	}

	// BWT Merging
	for(size_t i = 0; i < mergeVector.size(); i+=2)
	{
		if(i + 1 != mergeVector.size())
		{
			std::cout << "Merge1: " << mergeVector[i] << "\n";
			std::cout << "Merge2: " << mergeVector[i+1] << "\n";
		}
		else
		{
			std::cout << "Singleton: " << mergeVector[i] << "\n";
		}
	}

	return NULL;
}
