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
#include "BWTWriter.h"
#include "BWTReader.h"

#define WORKING_FILENAME "partialbwt."

// Definitions and structures
typedef uint16_t GAP_TYPE;

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

// Function declaration
MergeItem mergeBWTDisk(SeqReader* pReader, const MergeItem& item1, 
                       const MergeItem& item2, const std::string& outfile);

//
void writeMergedBWT(const std::string& outfile, const std::string& external_file, 
                    BWT* pBWTInternal, std::vector<GAP_TYPE>& gap_array);

//
void updateGapArray(const DNAString& w, const BWT* pBWTInternal, 
                    std::vector<GAP_TYPE>& gap_array);

//
inline void incrementGapArray(int64_t rank, std::vector<GAP_TYPE>& gap_array);

// The algorithm is as follows. We create M BWTs for subsets of 
// the input reads. These are created independently and written
// to disk. They are then merged either sequentially or pairwise
// to create the final BWT
void buildBWTDisk(const std::string& filename, const std::string& bwt_extension)
{
	size_t MAX_READS_PER_GROUP = 1;

	SeqReader* pReader = new SeqReader(filename);
	SeqRecord record;

	ReadTable* pCurrRT = new ReadTable;

	int groupID = 0;
	size_t numReadTotal = 0;

	MergeVector mergeVector;
	MergeItem mergeItem;
	mergeItem.start_index = 0;

	// Phase 1: Compute the initial BWTs
	bool done = false;
	while(!done)
	{
		done = !pReader->get(record);

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

				// Push the merge info
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
	delete pReader;

	// Phase 2: Pairwise merge the BWTs
	pReader = new SeqReader(filename);
	for(size_t i = 0; i < mergeVector.size(); i+=2)
	{
		if(i + 1 != mergeVector.size())
		{
			std::stringstream ss;
			ss << WORKING_FILENAME << groupID++ << bwt_extension;
			std::string temp_filename = ss.str();
			MergeItem merged = mergeBWTDisk(pReader, mergeVector[i], mergeVector[i+1], temp_filename);
		}
		else
		{
			std::cout << "Singleton: " << mergeVector[i] << "\n";
		}
	}
	delete pReader;
}

// Merge a pair of BWTs using disk storage
// Precondition: pReader is positioned at the start
// of the read block for item1
MergeItem mergeBWTDisk(SeqReader* pReader, const MergeItem& item1, 
                       const MergeItem& item2, const std::string& outfile)
{
	std::cout << "Merge1: " << item1 << "\n";
	std::cout << "Merge2: " << item2 << "\n";
	assert(item2.start_index == item1.end_index + 1);

	// Load the bwt of item2 into memory as the internal bwt
	BWT* pBWTInternal = new BWT(item2.filename);

	// Create the gap array, one value for each position in the internal bwt
	size_t gap_array_size = pBWTInternal->getBWLen();
	std::vector<GAP_TYPE> gap_array(gap_array_size, 0);

	// Calculate the rank of every read from item1.start_index to item1.end_index
	// and increment the gap counts
	size_t curr_idx = item1.start_index;
	SeqRecord record;
	while(curr_idx <= item1.end_index)
	{
		bool eof = !pReader->get(record);
		assert(!eof);
		DNAString& seq = record.seq;

		// Compute the ranks of all suffixes of seq
		updateGapArray(seq, pBWTInternal, gap_array);
		++curr_idx;
	}

	// Write the merged BWT to disk
	writeMergedBWT(outfile, item1.filename, pBWTInternal, gap_array);
	delete pBWTInternal;

	// Move the file pointer to the end of item2's reads.
	// This would be best done with a seek() but gzstream 
	// does not support this so we inefficiently just read through
	// the file. TODO: Change this.
	WARN_ONCE("Replace read through with seek()");
	while(curr_idx <= item2.end_index)
	{
		bool eof = !pReader->get(record);
		assert(!eof);
		++curr_idx;
	}
	assert(curr_idx = item2.end_index + 1);

	// Compute the merged item
	MergeItem merged;
	merged.start_index = item1.start_index;
	merged.end_index = item2.end_index;
	merged.filename = outfile;
	return merged;
}

// Write the header and BWStr of the merged product
void writeMergedBWT(const std::string& outfile, const std::string& external_file, 
                    BWT* pBWTInternal, std::vector<GAP_TYPE>& gap_array)
{
	BWTWriter* pWriter = new BWTWriter(outfile);
	BWTReader* pBWTDisk = new BWTReader(external_file);

	// Calculate and write header values
	size_t disk_strings;
	size_t disk_symbols;
	pBWTDisk->readHeader(disk_strings, disk_symbols);

	size_t total_strings = disk_strings + pBWTInternal->getNumStrings();
	size_t total_symbols = disk_symbols + pBWTInternal->getBWLen();
	pWriter->writeHeader(total_strings, total_symbols);
	
	// Calculate and write the actual string
	// The semantics of the gap array are that we need to write gap_array[i]
	// symbols to the stream before writing bwtInternal[i]
	size_t num_wrote = 0;
	for(size_t i = 0; i < pBWTInternal->getBWLen(); ++i)
	{
		size_t v = gap_array[i];
		std::cout << "Gap val: " << v << "\n";
		for(size_t j = 0; j < v; ++j)
		{
			char b = pBWTDisk->readBWChar();
			std::cout << "Read from disk: " << b << "\n";
			assert(b != '\n');
			pWriter->writeBWChar(b);
			++num_wrote;
		}

		pWriter->writeBWChar(pBWTInternal->getChar(i));
		++num_wrote;
	}
	
	assert(num_wrote == total_symbols);

	// Ensure we read the entire bw string from disk
	char last = pBWTDisk->readBWChar();
	assert(last == '\n');
	// Write a newline to finish the bwstr section
	pWriter->writeBWChar('\n');

	delete pWriter;
	delete pBWTDisk;
}


// Increment the gap array for each suffix of seq
void updateGapArray(const DNAString& w, const BWT* pBWTInternal, std::vector<GAP_TYPE>& gap_array)
{
	size_t l = w.length();
	int i = l - 1;

	// Compute the rank of the last character of seq. We consider $ to be implicitly
	// terminated by a $ character. The first rank calculated is for this and it is given
	// by the C(a) array in BWTInternal
	int64_t rank = pBWTInternal->getPC('$'); // always zero
	incrementGapArray(rank, gap_array);

	// Compute the starting rank for the last symbol of w
	char c = w.get(i);
	rank = pBWTInternal->getPC(c);
	incrementGapArray(rank, gap_array);
	printf("rank %c %d\n", '$', (int)rank);
	--i;

	// Iteratively compute the rest
	while(i >= 0)
	{
		char c = w.get(i);
		rank = pBWTInternal->getPC(c) + pBWTInternal->getOcc(c, rank - 1);
		incrementGapArray(rank, gap_array);
		printf("rank %c %d\n", c, (int)rank);
		--i;
	}
}

void incrementGapArray(int64_t rank, std::vector<GAP_TYPE>& gap_array)
{
	static size_t max_gap_count = std::numeric_limits<GAP_TYPE>::max();
	assert(gap_array[rank] < max_gap_count);
	++gap_array[rank];
}

