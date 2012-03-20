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

    // Start the haplotype generation process by finding reads containing the initial kmer
    HaplotypeReadVector positioned_reads;
    std::vector<std::string> initial_query_kmers;
    initial_query_kmers.push_back(m_initial_kmer_string);
    addPositionedReadsForKmers(m_initial_kmer_string, initial_query_kmers, &positioned_reads);

    int round = 0;
    bool extended = false;
    std::string consensus = "";
    do {
        extended = false;
        MultipleAlignment multiple_alignment = buildMultipleAlignment(positioned_reads);
        consensus = getConsensus(&multiple_alignment, 500);
        
        if(round % 10 == 0) 
        {
            printf("Round %d\n", round++);
            multiple_alignment.print(200);
        }

        // Attempt to extend the haplotype using newly found variant kmers if a consensus was computed
        if(!consensus.empty())
        {
            // Renumber the positioned reads to start from zero
            std::sort(positioned_reads.begin(), positioned_reads.end());
            
            int leftmost = positioned_reads.front().position;
            for(size_t i = 0; i < positioned_reads.size(); ++i)
                positioned_reads[i].position -= leftmost;

            std::vector<std::string> extension_kmer_vector = getExtensionKmers(consensus);
            if(!extension_kmer_vector.empty())
            {
                // Find the position of the extension kmer in the consensus
                addPositionedReadsForKmers(consensus, extension_kmer_vector, &positioned_reads);
                extended = true;
            }
        }
    } while(extended);

    if(!consensus.empty() && consensus.size() >= 105) 
    {
        std::cout << "Consensus: " << consensus << "\n";
        out_haplotypes.push_back(consensus);
    }

    return HBRC_OK;
}

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
std::vector<std::string> ReadCoherentHaplotypeBuilder::getExtensionKmers(const std::string& sequence)
{
    assert(sequence.size() >= m_parameters.kmer);
    std::vector<std::string> out_vector;
    for(size_t i = 0; i < sequence.size() - m_parameters.kmer + 1; ++i)
    {
        std::string kmer_sequence = sequence.substr(i, m_parameters.kmer);

        // This kmer has been used previously
        if(m_used_kmers.find(kmer_sequence) != m_used_kmers.end())
            continue;

        size_t var_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer_sequence, 
                                                                            m_parameters.pVariantBWT, 
                                                                            m_parameters.pVariantBWTCache);

        if(var_count < m_parameters.kmerThreshold)
            continue;

        // Check if this k-mer is present in the other base index
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer_sequence, 
                                                                             m_parameters.pBaseBWT, 
                                                                             m_parameters.pBaseBWTCache);

        if(base_count == 0)
            out_vector.push_back(kmer_sequence);
    }

    return out_vector;
}
