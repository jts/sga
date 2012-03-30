///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ReadCoherentHaplotypeBuilder - Construct candidate
// haplotypes from a pair of k-mer seeds.
//
#include "ReadCoherentHaplotypeBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "Profiler.h"
#include "HapgenUtil.h"
#include "multiple_alignment.h"

//
//
//
ReadCoherentHaplotypeBuilder::ReadCoherentHaplotypeBuilder()
{
}

//
ReadCoherentHaplotypeBuilder::~ReadCoherentHaplotypeBuilder()
{
}

// The source string is the string the bubble starts from
void ReadCoherentHaplotypeBuilder::setInitialHaplotype(const std::string& sequence)
{
    printf("\n\n***** Starting new haplotype %s\n", sequence.c_str());
    m_initial_kmer_string = sequence;
}

//
void ReadCoherentHaplotypeBuilder::setParameters(const GraphCompareParameters& parameters)
{
    m_parameters = parameters;
}

// Functor to sort the reads by the kmer position
struct KmerPositionSorter 
{
    KmerPositionSorter(const std::vector<size_t>& positionVector) : m_positions(positionVector) {}
    bool operator()(size_t i, size_t j) { return m_positions[i] > m_positions[j]; }
    const std::vector<size_t>& m_positions;
};

// Run the bubble construction process
HaplotypeBuilderReturnCode ReadCoherentHaplotypeBuilder::run(StringVector& out_haplotypes)
{
    PROFILE_FUNC("ReadCoherentHaplotypeBuilder::run")
    assert(m_parameters.pVariantBWT != NULL);
    assert(!m_initial_kmer_string.empty());

    //
    // Step 1. Generate an initial seed haplotype by aligning all the reads containing
    // the seed kmer against each other
    //
    std::string haplotype = generateSeedHaplotype(m_initial_kmer_string);
    if(haplotype.empty())
        return HBRC_OK;

    printf("Initial seed haplotype: %s\n", haplotype.c_str());

    // We expect the haplotype to eventually meet the sequence
    // that is in the base reads/reference
    bool right_rejoined = false;
    bool left_rejoined = false;
    int MAX_ROUNDS = 10;
    int round = 0;
    //
    // Right extension
    //
    bool continue_extension = true;
    while(continue_extension)
    {
        size_t length_before = haplotype.size();
        bool is_unambiguous = extendHaplotype(&haplotype);
        size_t length_after = haplotype.size();
        bool was_extended = length_after > length_before;
        
        // Check if haplotype rejoins base/reference sequence
        right_rejoined = checkAndRejoin(&haplotype, 10);
        continue_extension = round < MAX_ROUNDS && was_extended && is_unambiguous && !right_rejoined;
    }

    //
    // Left extension via reversing the haplotype
    //
    std::cout << "LEFT\n";
    haplotype = reverseComplement(haplotype);
    round = 0;
    continue_extension = true;
    while(continue_extension)
    {
        size_t length_before = haplotype.size();
        bool is_unambiguous = extendHaplotype(&haplotype);
        size_t length_after = haplotype.size();
        bool was_extended = length_after > length_before;
        
        // Check if haplotype rejoins base/reference sequence
        left_rejoined = checkAndRejoin(&haplotype, 10);
        continue_extension = round < MAX_ROUNDS && was_extended && is_unambiguous && !left_rejoined;
    }
    
    haplotype = reverseComplement(haplotype);
    printf("FINAL: %s\n", haplotype.c_str());

    if(haplotype.size() > 100 && left_rejoined && right_rejoined)
        out_haplotypes.push_back(haplotype);
    return HBRC_OK;
}

//
std::string ReadCoherentHaplotypeBuilder::generateSeedHaplotype(const std::string& kmer)
{
    // Extract all the reads that contain this kmer and order them by the position
    // of the kmer within the sequence
    HaplotypeReadVector positioned_reads;
    std::vector<std::string> initial_query_kmers;
    initial_query_kmers.push_back(kmer);
    addPositionedReadsForKmers(kmer, initial_query_kmers, &positioned_reads);
    assert(!positioned_reads.empty());

    // Compute a multiple alignment from the order of reads
    MultipleAlignment multiple_alignment = buildMultipleAlignmentUngapped(positioned_reads);
    multiple_alignment.print(200);

    // Compute an initial consensus for the multiple alignment
    // We limit the consensus to a window where all columns within the window have at least N depth
    std::vector<int> depth_vector;
    std::string initial = "";
    size_t total_columns = multiple_alignment.getNumColumns();
    for(size_t i = 0; i < total_columns; ++i)
    {
        SymbolCountVector symbol_counts = multiple_alignment.getSymbolCountVector(i);
        assert(!symbol_counts.empty());
        std::sort(symbol_counts.begin(), symbol_counts.end(), SymbolCount::countOrderDescending);
        
        // Calculate depth at this column
        int depth = 0;
        char max_base = '\0';
        for(size_t j = 0; j < symbol_counts.size(); ++j)
        {
            depth += symbol_counts[j].count;
            if(max_base == '\0')
                max_base = symbol_counts[j].symbol;
        }
        depth_vector.push_back(depth);
        initial.push_back(max_base);
    }

    // Calculate the maximal window in which the coverage is deep enough to make a call
    int MIN_DEPTH = m_parameters.kmerThreshold;
    std::vector<size_t> window_starts;
    std::vector<size_t> window_ends;

    bool in_window = false;
    for(size_t i = 0; i < depth_vector.size(); ++i)
    {
        if(!in_window && depth_vector[i] >= MIN_DEPTH)
        {
            window_starts.push_back(i);
            in_window = true;
        }

        if(in_window && depth_vector[i] < MIN_DEPTH)
        {
            window_ends.push_back(i);
            in_window = false;
        }
    }
    
    // end the last window
    if(in_window)
        window_ends.push_back(depth_vector.size());
    assert(window_starts.size() == window_ends.size());

    // Choose longest window
    size_t best_window_idx = 0;
    size_t best_window_size = 0;
    for(size_t i = 0; i < window_starts.size(); ++i)
    {
        size_t s = window_ends[i] - window_starts[i];
        if(s > best_window_size)
        {
            best_window_size = s;
            best_window_idx = i;
        }
    }

    if(best_window_size > 0)
    {
        printf("Best window [%zu %zu]\n", window_starts[best_window_idx], window_ends[best_window_idx]);
        std::string w = initial.substr(window_starts[best_window_idx], best_window_size);
        
        // Remove gaps
        std::string consensus;
        for(size_t i = 0; i < w.size(); ++i)
        {
            if(w[i] != '-')
                consensus.push_back(w[i]);
            else
                assert(false);
        }
        return w;
    }
    else
    {
        return "";
    }
}

static bool sortPairSecond(const std::pair<int, int>& a, const std::pair<int, int>& b)
{
    return a.second < b.second;
}

//
bool ReadCoherentHaplotypeBuilder::extendHaplotype(std::string* haplotype)
{
    // Recruit new reads that may overlap this haplotype
    std::set<std::string> kmer_set;
    size_t k = m_parameters.kmer;
    size_t nk = haplotype->size() - k + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = haplotype->substr(i, k);
        if(m_used_kmers.find(kmer) == m_used_kmers.end())
            kmer_set.insert(kmer);
    }

    std::vector<std::string> query_kmers(kmer_set.begin(), kmer_set.end());
    StringVector new_reads;
    getReadsForKmers(query_kmers, &new_reads);
    printf("Found %zu new reads from %zu kmers\n", new_reads.size(), query_kmers.size());

    // Filter out reads using following conditions:
    //   a) discard reads that do not have a proper suffix overlap of M bases with the haplotype
    //   b) discard reads that have more than N mismatches to the haplotype
    //   c) discard reads that do not match the haplotype at conflicting columns
    size_t CALL_DEPTH = m_parameters.kmerThreshold;
    int MAX_EDITS = 4;
    int MIN_OVERLAP = m_parameters.kmer;
    
    //
    std::vector<std::pair<int,int> > order_read_vector;

    // Overlap the new reads with the haplotype
    for(size_t i = 0; i < new_reads.size(); ++i)
    {
        SequenceOverlap overlap = Overlapper::computeOverlap(*haplotype, new_reads[i]);
        
        bool accept_read = true;
        // Check for proper overlap
        if(overlap.match[1].start != 0 || overlap.match[0].end != (int)haplotype->size() - 1)
            accept_read = false;
        
        // Check for overlap length
        if(overlap.getOverlapLength() < MIN_OVERLAP)
            accept_read = false;

        // Check for edits
        if(overlap.edit_distance > MAX_EDITS)
            accept_read = false;

        if(accept_read)
            order_read_vector.push_back(std::make_pair(i, overlap.match[0].start));
    }

    // Sort the ordering vector by the read alignment position along the haplotype
    std::sort(order_read_vector.begin(), order_read_vector.end(), sortPairSecond);

    // Build the multiple alignment by aligning each read with the previous read added
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("haplotype", *haplotype, "");
    std::string* previous = haplotype;
    for(size_t i = 0; i < order_read_vector.size(); ++i)
    {
        size_t idx = order_read_vector[i].first;

        SequenceOverlap overlap = Overlapper::computeOverlap(*previous, new_reads[idx], ungapped_params);

        // Ensure the overlap is a proper suffix overlap
        bool is_proper_0 = overlap.match[0].start > 0 && overlap.match[0].end == (int)previous->size() - 1;
        bool is_proper_1 = overlap.match[1].start == 0 && overlap.match[1].end < (int)previous->size() - 1;
        //bool has_gap = overlap.cigar.find_first_of("ID") != std::string::npos;
        if(is_proper_0 && is_proper_1)
        {
            //overlap.printAlignment(*previous, new_reads[idx]);
            multiple_alignment.addExtension("ext", new_reads[idx], "", overlap);
            previous = &new_reads[idx];
        }
    }

    // Filter the multiple alignment to achieve condition c
    multiple_alignment.print(200);
    multiple_alignment.filterByCount(CALL_DEPTH);
    if(!multiple_alignment.isValid())
        return false;
    multiple_alignment.print(200);

    // Compute the consensus sequence of the bases that extend the haplotype.
    // We only extend while there is an unambiguous extension and the sequence
    // is deep enough.
    
    // Find the first column after the last base of the haplotype
    size_t first_extend_column = 0;

    // Make sure there are no left-overhangs
    size_t total_columns = multiple_alignment.getNumColumns();
    assert(multiple_alignment.getSymbol(0, 0) != '\0');
    while(first_extend_column < total_columns && multiple_alignment.getSymbol(0, first_extend_column) != '\0')
        ++first_extend_column;
    assert(first_extend_column == total_columns || multiple_alignment.getSymbol(0, first_extend_column) == '\0');
    
    std::string extension = "";

    bool is_unambiguous = true;

    for(size_t i = first_extend_column; i < total_columns; ++i)
    {
        SymbolCountVector symbol_counts = multiple_alignment.getSymbolCountVector(i);
        assert(!symbol_counts.empty());
        std::sort(symbol_counts.begin(), symbol_counts.end(), SymbolCount::countOrderDescending);
        
        // Calculate depth at this column
        size_t total_depth = 0;
        char max_base = '\0';
        size_t max_depth = 0;
        size_t second_depth = 0;
        for(size_t j = 0; j < symbol_counts.size(); ++j)
        {
            total_depth += symbol_counts[j].count;

            // Record the first and second most frequent bases in this column
            if(max_depth == 0) 
            {
                max_depth = symbol_counts[j].count;
                max_base = symbol_counts[j].symbol;
            }
            else if(second_depth == 0)
            {
                second_depth = symbol_counts[j].count;
            }
        }

        // If the second-deepest base is greater than the MIN_CALL this 
        // column is ambiguous. Stop extension.
        if(second_depth >= CALL_DEPTH)
        {
            is_unambiguous = false;
            break;
        }

        if(total_depth >= CALL_DEPTH)
        {
            if(max_base != '-')
                extension.push_back(max_base);
        }
        else
        {
            break; // not enough depth to extend
        }
    }
    
    haplotype->append(extension);
    return is_unambiguous;
}

//
bool ReadCoherentHaplotypeBuilder::checkAndRejoin(std::string* haplotype, int max_distance)
{
    size_t k = m_parameters.kmer;
    if(haplotype->size() < k)
        return false;

    int nk = haplotype->size() - k + 1;
    int ki = nk - 1;
    while(ki >= 0 && ki >= nk - max_distance)
    {
        std::string ks = haplotype->substr(ki, k);
        // Check if present in base index
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(ks, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);
        if(base_count > 0)
        {
            // join found
            // trim haplotype to this kmer
            *haplotype = haplotype->substr(0, ki + k);
            return true;
        }
        ki -= 1;
    }

    // no join
    return false;
}

//
void ReadCoherentHaplotypeBuilder::addPositionedReadsForKmers(const std::string& consensus,
                                                              const std::vector<std::string>& kmer_vector,
                                                              HaplotypeReadVector* positioned_reads)
{
    SeqItemVector reads;
    SeqItemVector rc_reads;

    // Forward reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, m_parameters.kmer, false, 100000, &reads, NULL);

    // Reverse reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, m_parameters.kmer, true, 100000, &rc_reads, NULL);

    // Copy reads into the positioned read vector, initially with unset positions
    for(size_t i = 0; i < reads.size(); ++i)
    {
        std::string read_sequence = reads[i].seq.toString();

        // Find the position in the read of one of the query kmers
        size_t read_position = std::string::npos;
        size_t consensus_position = std::string::npos;
        std::vector<std::string>::const_iterator iter = kmer_vector.begin();
        while(read_position == std::string::npos && iter != kmer_vector.end())
        {
            read_position = read_sequence.find(*iter);
            consensus_position = consensus.find(*iter);
            ++iter;
        }

        assert(read_position != std::string::npos);
        assert(consensus_position != std::string::npos);

        HaplotypePositionedRead hpr = { reads[i].id, read_sequence, consensus_position - read_position };
        positioned_reads->push_back(hpr);
    }

    // reversed reads
    for(size_t i = 0; i < rc_reads.size(); ++i) 
    {
        std::string read_sequence = reverseComplement(rc_reads[i].seq.toString());

        // Find the position in the read of one of the query kmers
        size_t consensus_position = std::string::npos;
        size_t read_position = std::string::npos;
        std::vector<std::string>::const_iterator iter = kmer_vector.begin();
        while(read_position == std::string::npos && iter != kmer_vector.end())
        {
            read_position = read_sequence.find(*iter);
            consensus_position = consensus.find(*iter);
            ++iter;
        }

        assert(read_position != std::string::npos);
        assert(consensus_position != std::string::npos);

        HaplotypePositionedRead hpr = { rc_reads[i].id, read_sequence, consensus_position - read_position };
        positioned_reads->push_back(hpr);
    }

    // Sort by id and remove duplicates
    std::sort(positioned_reads->begin(), positioned_reads->end(), HaplotypePositionedRead::sortByID);
    HaplotypeReadVector::iterator iterator = std::unique(positioned_reads->begin(), positioned_reads->end(), HaplotypePositionedRead::equalByID);
    positioned_reads->resize(iterator - positioned_reads->begin());
    
    // record that we have used these kmers to extract reads.
    for(size_t i = 0; i < kmer_vector.size(); ++i)
        m_used_kmers.insert(kmer_vector[i]);
}

void ReadCoherentHaplotypeBuilder::getReadsForKmers(const std::vector<std::string>& kmer_vector,
                                                    std::vector<std::string>* out_reads)
{
    SeqItemVector reads;
    SeqItemVector rc_reads;
    size_t extraction_kmer = 31;

    // Forward reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, extraction_kmer, false, 100000, &reads, NULL);

    // Reverse reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, extraction_kmer, true, 100000, &rc_reads, NULL);

    // Copy reads into the positioned read vector, initially with unset positions
    for(size_t i = 0; i < reads.size(); ++i)
    {
        std::string read_sequence = reads[i].seq.toString();
        out_reads->push_back(read_sequence);
    }

    // reversed reads
    for(size_t i = 0; i < rc_reads.size(); ++i) 
    {
        std::string read_sequence = reverseComplement(rc_reads[i].seq.toString());
        out_reads->push_back(read_sequence);
    }

    // record that we have used these kmers to extract reads.
    for(size_t i = 0; i < kmer_vector.size(); ++i)
        m_used_kmers.insert(kmer_vector[i]);
}

//
MultipleAlignment ReadCoherentHaplotypeBuilder::buildMultipleAlignmentUngapped(HaplotypeReadVector& positioned_reads) const
{
    // Sort the haplotype reads into ascending order by position
    std::sort(positioned_reads.begin(), positioned_reads.end());
    assert(!positioned_reads.empty());

    MultipleAlignment multiple_alignment;

    // Add leftmost sequence as the first read of the multiple alignment
    multiple_alignment.addBaseSequence("base", positioned_reads[0].sequence, "");
    
    // Add remaining sequences
    for(size_t i = 1; i < positioned_reads.size(); ++i)
    {
        // Craft an overlap object for this sequence, relative to the last sequence added
        size_t prev = i - 1;

        size_t p0 = positioned_reads[prev].position;
        size_t p1 = positioned_reads[i].position;

        SequenceOverlap overlap;
        overlap.match[0].start = p1 - p0;
        overlap.match[0].end = positioned_reads[prev].sequence.length() - 1;

        overlap.match[1].start = 0;
        overlap.match[1].end = positioned_reads[i].sequence.length() - 1 - (p1 - p0);

        assert(overlap.match[0].length() == overlap.match[1].length());
        size_t overlap_length = overlap.match[0].length();
        std::stringstream cigar_builder;
        cigar_builder << overlap_length << "M";
        overlap.cigar = cigar_builder.str();
        overlap.score = 2*overlap_length;
        overlap.edit_distance = -1;
        overlap.total_columns = overlap_length;

        multiple_alignment.addExtension("seq", positioned_reads[i].sequence, "", overlap);
    }

    return multiple_alignment;
}
