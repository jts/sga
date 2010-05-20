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
#include "SAWriter.h"
#include "SAReader.h"

// Definitions and structures
typedef uint32_t GAP_TYPE;
static const bool USE_GZ = false;
struct MergeItem
{
    int64_t start_index;
    int64_t end_index;
    std::string reads_filename;
    std::string bwt_filename;
    std::string sai_filename;

    friend std::ostream& operator<<(std::ostream& out, const MergeItem& item)
    {
        out << "[" << item.start_index << "," << item.end_index << "] " << item.bwt_filename;
        out << " " << item.sai_filename << " " << item.reads_filename;
        return out;
    }
};
typedef std::vector<MergeItem> MergeVector;

// Function declarations
int64_t merge(SeqReader* pReader, 
              const MergeItem& item1, const MergeItem& item2, 
              const std::string& bwt_outname, const std::string& sai_outname,
              bool doReverse);


//
void writeMergedIndex(const BWT* pBWTInternal, const MergeItem& externalItem, 
                      const MergeItem& internalItem, const std::string& bwt_outname,
                      const std::string& sai_outname, const std::vector<GAP_TYPE>& gap_array);

//
void updateGapArray(const DNAString& w, const BWT* pBWTInternal, 
                    std::vector<GAP_TYPE>& gap_array);

//
inline void incrementGapArray(int64_t rank, std::vector<GAP_TYPE>& gap_array);
std::string makeTempName(const std::string& prefix, int id, const std::string& extension);
std::string makeFilename(const std::string& prefix, const std::string& extension);


// The algorithm is as follows. We create M BWTs for subsets of 
// the input reads. These are created independently and written
// to disk. They are then merged either sequentially or pairwise
// to create the final BWT
void buildBWTDisk(const std::string& in_filename, const std::string& out_prefix, 
                  const std::string& bwt_extension, const std::string& sai_extension,
                  bool doReverse)
{
    size_t MAX_READS_PER_GROUP = 200000;

    SeqReader* pReader = new SeqReader(in_filename);
    SeqRecord record;

    int groupID = 0;
    size_t numReadTotal = 0;

    MergeVector mergeVector;
    MergeItem mergeItem;
    mergeItem.start_index = 0;

    // Phase 1: Compute the initial BWTs
    ReadTable* pCurrRT = new ReadTable;
    bool done = false;
    while(!done)
    {
        done = !pReader->get(record);

        if(!done)
        {
            // the read is valid
            SeqItem item = record.toSeqItem();
            if(doReverse)
                item.seq.reverse();
            pCurrRT->addRead(item);
            ++numReadTotal;
        }

        if(pCurrRT->getCount() >= MAX_READS_PER_GROUP || (done && pCurrRT->getCount() > 0))
        {
            // Compute the SA and BWT for this group
            SuffixArray* pSA = new SuffixArray(pCurrRT);
            BWT* pBWT = new BWT(pSA, pCurrRT);

            // Write the BWT to disk                
            std::string bwt_temp_filename = makeTempName(out_prefix, groupID, bwt_extension);
            pBWT->write(bwt_temp_filename);

            std::string sai_temp_filename = makeTempName(out_prefix, groupID, sai_extension);
            pSA->writeIndex(sai_temp_filename);

            // Push the merge info
            mergeItem.end_index = numReadTotal - 1;
            mergeItem.reads_filename = in_filename;
            mergeItem.bwt_filename = bwt_temp_filename;
            mergeItem.sai_filename = sai_temp_filename;
            mergeVector.push_back(mergeItem);

            // Cleanup
            delete pBWT;
            delete pSA;

            // Start the new group
            mergeItem.start_index = numReadTotal;
            ++groupID;
            pCurrRT->clear();
        }
    }
    delete pCurrRT;
    delete pReader;

    // Phase 2: Pairwise merge the BWTs
    int round = 1;
    MergeVector nextMergeRound;
    while(mergeVector.size() > 1)
    {
        std::cout << "Starting round " << round << "\n";
        pReader = new SeqReader(in_filename);
        for(size_t i = 0; i < mergeVector.size(); i+=2)
        {
            if(i + 1 != mergeVector.size())
            {
                std::string bwt_merged_name = makeTempName(out_prefix, groupID, bwt_extension);
                std::string sai_merged_name = makeTempName(out_prefix, groupID, sai_extension);

                MergeItem item1 = mergeVector[i];
                MergeItem item2 = mergeVector[i+1];

                // Perform the actual merge
                int64_t curr_idx = merge(pReader, item1, item2, 
                                         bwt_merged_name, sai_merged_name, 
                                         doReverse);

                // pReader now points to the end of item1's block of 
                // reads. Skip item2's reads
                assert(curr_idx == item2.start_index);
                while(curr_idx <= item2.end_index)
                {
                    bool eof = !pReader->get(record);
                    assert(!eof);
                    ++curr_idx;
                }

                // Create the merged mergeItem to use in the next round
                MergeItem merged;
                merged.start_index = item1.start_index;
                merged.end_index = item2.end_index;
                merged.bwt_filename = bwt_merged_name;
                merged.sai_filename = sai_merged_name;
                nextMergeRound.push_back(merged);

                // Done with the temp files, remove them
                unlink(item1.bwt_filename.c_str());
                unlink(item2.bwt_filename.c_str());
                unlink(item1.sai_filename.c_str());
                unlink(item2.sai_filename.c_str());

                ++groupID;
            }
            else
            {
                // Singleton, pass through to the next round
                nextMergeRound.push_back(mergeVector[i]);
            }
        }
        delete pReader;
        mergeVector.clear();
        mergeVector.swap(nextMergeRound);
        ++round;
    }
    assert(mergeVector.size() == 1);

    // Done, rename the files to their final name
    std::stringstream bwt_ss;
    bwt_ss << out_prefix << bwt_extension << (USE_GZ ? ".gz" : "");
    std::string bwt_final_filename = bwt_ss.str();
    rename(mergeVector.front().bwt_filename.c_str(), bwt_final_filename.c_str());

    std::stringstream sai_ss;
    sai_ss << out_prefix << sai_extension << (USE_GZ ? ".gz" : "");
    std::string sai_final_filename = sai_ss.str();
    rename(mergeVector.front().sai_filename.c_str(), sai_final_filename.c_str());
}

// Merge the indices for the two independent sets of reads in readsFile1 and readsFile2
void mergeIndependentIndices(const std::string& readsFile1, const std::string& readsFile2, 
                             const std::string& bwt_extension, const std::string& sai_extension, 
                             bool doReverse)
{
    MergeItem item1;
    std::string prefix1 = stripFilename(readsFile1);
    item1.reads_filename = readsFile1;
    item1.bwt_filename = makeFilename(prefix1, bwt_extension);
    item1.sai_filename = makeFilename(prefix1, sai_extension);
    item1.start_index = 0;
    item1.end_index = -1; // this tells merge to read the entire file

    MergeItem item2;
    std::string prefix2 = stripFilename(readsFile2);
    item2.reads_filename = readsFile2;
    item2.bwt_filename = makeFilename(prefix2, bwt_extension);
    item2.sai_filename = makeFilename(prefix2, sai_extension);
    item2.start_index = 0;
    item2.end_index = -1;

    // Prepare the input filehandle
    SeqReader* pReader = new SeqReader(item1.reads_filename);
    
    // Build the outnames
    std::string merged_prefix = prefix1 + "." + prefix2;
    std::string bwt_merged_name = makeFilename(merged_prefix, bwt_extension);
    std::string sai_merged_name = makeFilename(merged_prefix, sai_extension);

    // Perform the actual merge
    merge(pReader, item1, item2, bwt_merged_name, sai_merged_name, doReverse);
    delete pReader;
}

// Merge two readsFiles together
void mergeReadFiles(const std::string& readsFile1, const std::string& readsFile2, const std::string& outFile)
{
    // If the outfile is the empty string, append the reads in readsFile2 into readsFile1
    // otherwise cat the files together
    std::ostream* pWriter;
    if(outFile.empty())
    {
        pWriter = createWriter(readsFile1, std::ios_base::out | std::ios_base::app);
    }
    else
    {
        pWriter = createWriter(outFile);

        // Copy reads1 to the outfile
        SeqReader reader(readsFile1);
        SeqRecord record;
        while(reader.get(record))
            record.write(*pWriter);
    }

    // Copy reads2 to writer
    SeqReader reader(readsFile2);
    SeqRecord record;
    while(reader.get(record))
        record.write(*pWriter);
    delete pWriter;
}

// Merge a pair of BWTs using disk storage
// Precondition: pReader is positioned at the start of the read block for item1
int64_t merge(SeqReader* pReader,
              const MergeItem& item1, const MergeItem& item2,
              const std::string& bwt_outname, const std::string& sai_outname,
              bool doReverse)
{
    std::cout << "Merge1: " << item1 << "\n";
    std::cout << "Merge2: " << item2 << "\n";

    // Load the bwt of item2 into memory as the internal bwt
    BWT* pBWTInternal = new BWT(item2.bwt_filename);

    // Create the gap array
    size_t gap_array_size = pBWTInternal->getBWLen() + 1;
    std::vector<GAP_TYPE> gap_array(gap_array_size, 0);

    // Calculate the rank of every read from item1.start_index to item1.end_index
    // and increment the gap counts
    int64_t curr_idx = item1.start_index;
    SeqRecord record;
    while(pReader->get(record))
    {
        DNAString& seq = record.seq;
        if(doReverse)
            seq.reverse();

        // Compute the ranks of all suffixes of seq
        updateGapArray(seq, pBWTInternal, gap_array);
        ++curr_idx;

        // If the end index is -1 we want to process every read,
        // otherwise stop once we are outside this block's range
        if(item1.end_index != -1 && curr_idx > item1.end_index)
            break;
    }
    assert(item1.end_index == -1 || (curr_idx == item1.end_index + 1 && curr_idx == item2.start_index));

    // Write the merged BWT/SAI to disk
    writeMergedIndex(pBWTInternal, item1, item2, bwt_outname, sai_outname, gap_array);
    delete pBWTInternal;
    return curr_idx;
}

// Merge the internal and external BWTs and the SAIs
void writeMergedIndex(const BWT* pBWTInternal, const MergeItem& externalItem, 
                      const MergeItem& internalItem, const std::string& bwt_outname,
                      const std::string& sai_outname, const std::vector<GAP_TYPE>& gap_array)
{
    BWTWriter bwtWriter(bwt_outname);
    BWTReader bwtExtReader(externalItem.bwt_filename);
    
    SAWriter saiWriter(sai_outname);
    SAReader saiExtReader(externalItem.sai_filename);
    SAReader saiIntReader(internalItem.sai_filename);

    // Calculate and write header values
    size_t disk_strings;
    size_t disk_symbols;
    BWFlag flag;
    bwtExtReader.readHeader(disk_strings, disk_symbols, flag);

    size_t total_strings = disk_strings + pBWTInternal->getNumStrings();
    size_t total_symbols = disk_symbols + pBWTInternal->getBWLen();
    bwtWriter.writeHeader(total_strings, total_symbols, BWF_NOFMI);
    
    // Discard the first two elements of each sai
    size_t discard1, discard2;
    saiExtReader.readHeader(discard1, discard2);
    saiIntReader.readHeader(discard1, discard2);

    // Write the header of the SAI which is just the number of strings and elements in the SAI
    saiWriter.writeHeader(total_strings, total_strings);

    // Calculate and write the actual string
    // The semantics of the gap array are that we need to write gap_array[i]
    // symbols to the stream before writing bwtInternal[i]
    // Each time a '$' symbol is read, it signals the end of some read. We
    // output one element of the sai from corresponding internal or external
    // sai file.
    size_t num_bwt_wrote = 0;
    size_t num_sai_wrote = 0;
    for(size_t i = 0; i < gap_array.size(); ++i)
    {
        size_t v = gap_array[i];
        for(size_t j = 0; j < v; ++j)
        {
            char b = bwtExtReader.readBWChar();
            assert(b != '\n');
            bwtWriter.writeBWChar(b);
            ++num_bwt_wrote;
            
            if(b == '$')
            {
                // The external indices are correct and only need to be copied
                SAElem e = saiExtReader.readElem(); 
                saiWriter.writeElem(e);
                ++num_sai_wrote;
            }
        }
        
        // If this is the last entry in the gap array, do not output a symbol from
        // the internal BWT
        if(i != pBWTInternal->getBWLen())
        {
            char b = pBWTInternal->getChar(i);
            bwtWriter.writeBWChar(b);
            ++num_bwt_wrote;

            if(b == '$')
            {
                // The internal indices need to be offset
                // by the number of strings in the external collection
                SAElem e = saiIntReader.readElem(); 

                uint64_t id = e.getID();
                id += disk_strings;
                e.setID(id);
                
                saiWriter.writeElem(e);
                ++num_sai_wrote;
            }
        }
    }
    assert(num_bwt_wrote == total_symbols);
    assert(num_sai_wrote == total_strings);

    // Ensure we read the entire bw string from disk
    char last = bwtExtReader.readBWChar();
    assert(last == '\n');
    // Write a newline to finish the bwstr section
    bwtWriter.writeBWChar('\n');
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
    --i;

    // Iteratively compute the remaining ranks
    while(i >= 0)
    {
        char c = w.get(i);
        rank = pBWTInternal->getPC(c) + pBWTInternal->getOcc(c, rank - 1);
        //std::cout << "c: " << c << " rank: " << rank << "\n";
        incrementGapArray(rank, gap_array);
        --i;
    }
}

//
void incrementGapArray(int64_t rank, std::vector<GAP_TYPE>& gap_array)
{
    static size_t max_gap_count = std::numeric_limits<GAP_TYPE>::max();
    assert(gap_array[rank] < max_gap_count);
    assert(rank < (int64_t)gap_array.size());
    ++gap_array[rank];
}

//
std::string makeTempName(const std::string& prefix, int id, const std::string& extension)
{
    std::stringstream ss;
    ss << prefix << ".temp-" << id << extension << (USE_GZ ? ".gz" : "");
    return ss.str();
}

//
std::string makeFilename(const std::string& prefix, const std::string& extension)
{
    std::stringstream ss;
    ss << prefix << extension << (USE_GZ ? ".gz" : "");
    return ss.str();
}
