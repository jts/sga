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
ReadCoherentHaplotypeBuilder::ReadCoherentHaplotypeBuilder() : m_pBWT(NULL), m_kmer(0), m_maxEditDistance(1)
{
}

//
ReadCoherentHaplotypeBuilder::~ReadCoherentHaplotypeBuilder()
{
}

//
void ReadCoherentHaplotypeBuilder::setKmer(size_t k)
{
    m_kmer = k;
    m_kmerThreshold = 1;
}

// The source string is the string the bubble starts from
void ReadCoherentHaplotypeBuilder::setInitialHaplotype(const std::string& sequence)
{
    printf("Starting new haplotype %s\n", sequence.c_str());
    m_initial_kmer_string = sequence;
}

// The source index is the index that the contains the source string
void ReadCoherentHaplotypeBuilder::setIndex(const BWT* pBWT, const BWTIntervalCache* pCache, 
                                            const SampledSuffixArray* pSSA)
{
    m_pBWT = pBWT;
    m_pIntervalCache = pCache;
    m_pSSA = pSSA;
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
    assert(m_pBWT != NULL);
    assert(!m_initial_kmer_string.empty());

    // Extract reads with the input k-mer
    std::vector<std::string> read_vector;
    getReadsWithKmer(m_initial_kmer_string, &read_vector);

    // Infer the position of the first base of each read relative to the underlying haplotype
    HaplotypeReadVector positioned_reads;
    for(size_t i = 0; i < read_vector.size(); ++i)
    {
        size_t pos = read_vector[i].find(m_initial_kmer_string);
        HaplotypePositionedRead hpr = { read_vector[i], -pos };
        positioned_reads.push_back(hpr);
    }

    MultipleAlignment multiple_alignment = buildMultipleAlignment(positioned_reads);
    multiple_alignment.print(200);
    
    std::string consensus = getConsensus(&multiple_alignment, 5);
    if(!consensus.empty() && consensus.size() >= 110) 
    {
        std::cout << "Consensus: " << consensus << "\n";
        out_haplotypes.push_back(consensus);
    }
#if 0
    // Sort reads by kmer position
    KmerPositionSorter sorter(read_kmer_positions);
    std::sort(read_indices.begin(), read_indices.end(), sorter);

    // Build a multiple alignment
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence("base", read_vector[read_indices.front()], "");
    size_t previous_added = 0;
    size_t total_added = 0;
    for(size_t i = 1; i < read_indices.size(); ++i)
    { 
        size_t read_idx_1 = read_indices[previous_added];
        size_t read_idx_2 = read_indices[i];

        SequenceOverlap overlap = Overlapper::extendMatch(read_vector[read_idx_1], read_vector[read_idx_2], read_kmer_positions[read_idx_1], read_kmer_positions[read_idx_2], 1);

        // Check if this is a valid overlap.
        bool is_proper_overlap = overlap.match[1].end < overlap.length[1] - 1 && overlap.match[1].start == 0;
        bool is_within_tolerance = overlap.edit_distance <= 1;
        bool is_gapless = overlap.total_columns == overlap.match[0].length() && overlap.total_columns == overlap.match[1].length();

        if(is_proper_overlap && is_within_tolerance && is_gapless)
        {
            std::cout << "ACCEPTED\n";
            overlap.printAlignment(read_vector[read_idx_1], read_vector[read_idx_2]);

            multiple_alignment.addExtension("test", read_vector[read_idx_2], "", overlap);
            previous_added = i;
            total_added += 1;
        }
        else
        {
            std::cout << "REJECTED\n";
            overlap.printAlignment(read_vector[read_idx_1], read_vector[read_idx_2]);
        }
    }

    printf("Successfully added %zu/%zu to multiple alignment\n", total_added, read_vector.size());
    multiple_alignment.print(200);

    // Compute consensus
    std::string consensus = "";
    if(total_added >= 3) 
    {
        consensus = getConsensus(&multiple_alignment);
        out_haplotypes.push_back(consensus);
    }
#endif
    return HBRC_OK;
}

//
MultipleAlignment ReadCoherentHaplotypeBuilder::buildMultipleAlignment(HaplotypeReadVector& positioned_reads) const
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


//
std::string ReadCoherentHaplotypeBuilder::getConsensus(MultipleAlignment* multiple_alignment, int max_differences) const
{
    std::string consensus = "";
    size_t total_columns = multiple_alignment->getNumColumns();
    for(size_t i = 0; i < total_columns; ++i)
    {
        SymbolCountVector symbol_counts = multiple_alignment->getSymbolCountVector(i);
        assert(!symbol_counts.empty());

        std::sort(symbol_counts.begin(), symbol_counts.end(), SymbolCount::countOrderDescending);
        // If this column has a maximum of symbols (2 or more bases with count >= max_differences)
        // then we saw this is a divergent column and return no haplotype.
        if(symbol_counts.size() >= 2)
        {
            if(symbol_counts[1].count >= max_differences)
                return "";
        }

        // Not divergent, add the most frequent base to the consensus
        size_t idx = 0;
        while(symbol_counts[idx].symbol == '-' || symbol_counts[idx].symbol == '\0')
            idx++;
        consensus.append(1, symbol_counts[idx].symbol);
    }
    return consensus;
}

//
void ReadCoherentHaplotypeBuilder::getReadsWithKmer(const std::string& kmer, std::vector<std::string>* out_reads)
{
    SeqItemVector variant_reads;
    SeqItemVector variant_rc_reads;

    std::vector<std::string> lookup_vector;
    lookup_vector.push_back(kmer);

    // Forward reads
    HapgenUtil::extractHaplotypeReads(lookup_vector, m_pBWT, m_pIntervalCache, m_pSSA, 
                                      m_kmer, false, 100000, &variant_reads, NULL);

    // Reverse reads
    HapgenUtil::extractHaplotypeReads(lookup_vector, m_pBWT, m_pIntervalCache, m_pSSA, 
                                      m_kmer, true, 100000, &variant_rc_reads, NULL);

    // Move read sequences to the vector of reads
    for(size_t i = 0; i < variant_reads.size(); ++i)
        out_reads->push_back(variant_reads[i].seq.toString());

    for(size_t i = 0; i < variant_rc_reads.size(); ++i)
        out_reads->push_back(reverseComplement(variant_rc_reads[i].seq.toString()));
}
