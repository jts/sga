///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapHaplotypeBuilder - Construct candidate
// haplotypes from a pair of k-mer seeds.
//
#include "OverlapHaplotypeBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "Profiler.h"
#include "HapgenUtil.h"
#include "multiple_alignment.h"

//
//
//
OverlapHaplotypeBuilder::OverlapHaplotypeBuilder()
{
}

//
OverlapHaplotypeBuilder::~OverlapHaplotypeBuilder()
{
}

// The source string is the string the bubble starts from
void OverlapHaplotypeBuilder::setInitialHaplotype(const std::string& sequence)
{
    printf("\n\n***** Starting new haplotype %s\n", sequence.c_str());
    m_initial_kmer_string = sequence;
}

//
void OverlapHaplotypeBuilder::setParameters(const GraphCompareParameters& parameters)
{
    m_parameters = parameters;
}

// Run the bubble construction process
HaplotypeBuilderReturnCode OverlapHaplotypeBuilder::run(StringVector& out_haplotypes)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::run")
    assert(m_parameters.pVariantBWT != NULL);
    assert(!m_initial_kmer_string.empty());

    // Start the haplotype generation process by finding reads containing the initial kmer
    StringVector query_kmers(1, m_initial_kmer_string);
    StringVector reads;
    getReadsForKmers(query_kmers, &reads);
    
    printf("Starting processing with %zu reads\n", reads.size());
    StringVector ordered_reads;
    orderReadsInitial(m_initial_kmer_string, reads, &ordered_reads);
    MultipleAlignment initial_ma = buildMultipleAlignment(ordered_reads);
    initial_ma.print(200);

    // Extend
    std::string consensus = "";
    bool done = false;
    while(!done)
    {
        StringVector extension_kmers = getExtensionKmers(reads);
        StringVector incoming_reads;
        getReadsForKmers(extension_kmers, &incoming_reads);
        orderReadsExtended(incoming_reads, &ordered_reads);
        removeDuplicates(&ordered_reads);

        printf("Extended to %zu reads\n", incoming_reads.size());
        MultipleAlignment extended_ma = buildMultipleAlignment(ordered_reads);
        extended_ma.print(400);

        // Try to compute a consensus
        consensus = getConsensus(&extended_ma, 3, 1000);

        // Try to find left/right anchors on the consensus
        int clip_k = m_parameters.kmer;
        int num_kmers = consensus.size() - clip_k + 1;
        int left_clip = num_kmers;
        int right_clip = -1;

        for(int i = 0; i < num_kmers; ++i)
        {
            std::string kmer_sequence = consensus.substr(i, clip_k);

            // Check whether the kmer exists in the base reads/reference
            size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer_sequence,
                                                                      m_parameters.pBaseBWT, 
                                                                      m_parameters.pBaseBWTCache);
            if(base_count > 0 && i < left_clip)
                left_clip = i;
            if(base_count > 0 && i > right_clip)
                right_clip = i;
        }
        
        printf("Clip range: [%d %d]\n", left_clip, right_clip);
        if(left_clip == num_kmers || right_clip == -1)
        {
            consensus = "";
            done = false;
        }
        else
        {
            consensus = consensus.substr(left_clip, right_clip - left_clip + clip_k);
            done = true;
        }
    }

    if(!consensus.empty())
        out_haplotypes.push_back(consensus);

    return HBRC_OK;
}

//
void OverlapHaplotypeBuilder::getReadsForKmers(const StringVector& kmer_vector, StringVector* reads)
{
    SeqItemVector fwd_si;
    SeqItemVector rev_si;

    // Forward reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, m_parameters.kmer, false, 100000, &fwd_si, NULL);

    // Reverse reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, m_parameters.kmer, true, 100000, &rev_si, NULL);

    // Copy reads into the positioned read vector, initially with unset positions
    for(size_t i = 0; i < fwd_si.size(); ++i)
    {
        std::string read_sequence = fwd_si[i].seq.toString();
        reads->push_back(read_sequence);
    }

    // reversed reads
    for(size_t i = 0; i < rev_si.size(); ++i) 
    {
        std::string read_sequence = reverseComplement(rev_si[i].seq.toString());
        reads->push_back(read_sequence);
    }

    // record that we have used these kmers to extract reads.
    for(size_t i = 0; i < kmer_vector.size(); ++i)
        m_used_kmers.insert(kmer_vector[i]);
}

static bool sortPairSecondAscending(const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b)
{
    return a.second > b.second;
}

// Order the list of reads based on the fact that they all share the same kmer
void OverlapHaplotypeBuilder::orderReadsInitial(const std::string& initial_kmer, const StringVector& reads, StringVector* ordered_vector)
{
    std::vector<std::pair<size_t, size_t> > read_kmer_vector;;
    for(size_t i = 0; i < reads.size(); ++i)
    {
        size_t pos = reads[i].find(initial_kmer);
        assert(pos != std::string::npos);
        read_kmer_vector.push_back(std::make_pair(i, pos));
    }

    // Sort the vector by the second element
    std::sort(read_kmer_vector.begin(), read_kmer_vector.end(), sortPairSecondAscending);

    // Insert the reads into the ordered list
    assert(ordered_vector->empty());
    for(size_t i = 0; i < read_kmer_vector.size(); ++i)
        ordered_vector->push_back(reads[read_kmer_vector[i].first]);

    // Build multiple alignment
    assert(!ordered_vector->empty());
}

// Insert the incoming reads into the ordered read list
void OverlapHaplotypeBuilder::orderReadsExtended(const StringVector& incoming_reads, StringVector* ordered_vector)
{
    int MIN_OVERLAP = 30;
    // Turn the vector into a list for efficient insertion
    StringList ordered_list(ordered_vector->begin(), ordered_vector->end());

    // Find the insertion position for every read
    // It is a precondition that every read has a significant overlap
    // with a read that is already ordered
    for(size_t i = 0; i < incoming_reads.size(); ++i)
    {
        StringList::iterator insert_iterator = ordered_list.begin();
        while(insert_iterator != ordered_list.end())
        {
            SequenceOverlap overlap = Overlapper::computeOverlap(incoming_reads[i], *insert_iterator);
            bool is_left_overlap = overlap.match[0].start == 0;
            bool is_right_overlap = overlap.match[0].end == (int)incoming_reads[i].size() - 1;
            bool is_contain = !is_left_overlap && !is_right_overlap;
            assert(!is_contain);
            bool is_full_match = is_left_overlap && is_right_overlap;

            // Check if we can insert here
            if(is_full_match || (is_right_overlap && overlap.getOverlapLength() >= MIN_OVERLAP))
                break;

            insert_iterator++;
        }
        
        ordered_list.insert(insert_iterator, incoming_reads[i]);
    }
    
    // Turn the list back into a vector for output
    ordered_vector->clear();
    ordered_vector->insert(ordered_vector->end(), ordered_list.begin(), ordered_list.end());
}

//
static bool sortPairFirst(const std::pair<size_t, std::string>& a, const std::pair<size_t, std::string>& b)
{
    return a.first < b.first;
}

static bool sortPairSecond(const std::pair<size_t, std::string>& a, const std::pair<size_t, std::string>& b)
{
    return a.second < b.second;
}

static bool equalPairSecond(const std::pair<size_t, std::string>& a, const std::pair<size_t, std::string>& b)
{
    return a.second == b.second;
}

//
void OverlapHaplotypeBuilder::removeDuplicates(StringVector* ordered_vector)
{
    std::vector<std::pair<size_t, std::string> > index_sequence_vector;
    for(size_t i = 0; i < ordered_vector->size(); ++i)
        index_sequence_vector.push_back(std::make_pair(i, ordered_vector->at(i)));

    // Sort by sequence
    std::sort(index_sequence_vector.begin(), index_sequence_vector.end(), sortPairSecond);

    // Unique by sequence
    std::vector<std::pair<size_t, std::string> >::iterator new_end = std::unique(index_sequence_vector.begin(), index_sequence_vector.end(), equalPairSecond);
    index_sequence_vector.resize(new_end - index_sequence_vector.begin());

    // Sort by position
    std::sort(index_sequence_vector.begin(), index_sequence_vector.end(), sortPairFirst);

    // Copy back into the ordered vector
    ordered_vector->clear();
    for(size_t i = 0; i < index_sequence_vector.size(); ++i)
        ordered_vector->push_back(index_sequence_vector[i].second);
}

//
MultipleAlignment OverlapHaplotypeBuilder::buildMultipleAlignment(const StringVector& ordered_vector)
{
    MultipleAlignment multiple_alignment;

    // Insert first sequence
    StringVector::const_iterator prev_iterator = ordered_vector.begin();
    multiple_alignment.addBaseSequence("base", *prev_iterator, "");

    // Insert remaining sequences
    StringVector::const_iterator read_iterator = prev_iterator;
    read_iterator++;
    while(read_iterator != ordered_vector.end())
    {
        SequenceOverlap overlap = Overlapper::computeOverlap(*prev_iterator, *read_iterator, ungapped_params);
        overlap.printAlignment(*prev_iterator, *read_iterator);
        multiple_alignment.addExtension("seq", *read_iterator, "", overlap);
        prev_iterator = read_iterator;
        read_iterator++;
    }

    return multiple_alignment;
}

//
StringVector OverlapHaplotypeBuilder::getExtensionKmers(const StringVector& reads)
{
    StringVector out_kmers;
    size_t k = m_parameters.kmer;
    for(size_t i = 0; i < reads.size(); ++i)
    {
        size_t nk = reads[i].size() - k + 1;

        // Vector to track the kmers that are unique in this read
        StringVector read_unique_kmers;
        bool read_has_base_kmer = false;
        for(size_t j = 0; j < nk; ++j)
        {
            std::string kmer_sequence = reads[i].substr(j, k);

            // Check whether the kmer exists in the base reads/reference
            /*
            size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer_sequence, 
                                                                                 m_parameters.pBaseBWT, 
                                                                                 m_parameters.pBaseBWTCache);
        
            if(base_count > 0)
            {
                read_has_base_kmer = true;
                break;
            }
            */
            // Skip kmers that have been used before
            bool used = m_used_kmers.find(kmer_sequence) != m_used_kmers.end();
            if(used)
                continue;
        
           // Check that this kmer has occured enough times in the variant reads
           size_t var_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer_sequence, 
                                                                               m_parameters.pVariantBWT, 
                                                                               m_parameters.pVariantBWTCache);

           //
           if(var_count >= 1)
               read_unique_kmers.push_back(kmer_sequence);

        }
        
        // Do not extend using reads that rejoin the base graph
        if(read_has_base_kmer)
            read_unique_kmers.clear();
        else
            out_kmers.insert(out_kmers.end(), read_unique_kmers.begin(), read_unique_kmers.end());
    }

    // Remove kmers that are present in multiple reads
    std::sort(out_kmers.begin(), out_kmers.end());
    StringVector::iterator new_end = std::unique(out_kmers.begin(), out_kmers.end());
    out_kmers.resize(new_end - out_kmers.begin());
    return out_kmers;
}

//
std::string OverlapHaplotypeBuilder::getConsensus(MultipleAlignment* multiple_alignment, int min_call_coverage, int max_differences) const
{
    std::string initial_consensus = "";
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
        
        // Insert no-call if the consensus base does not have enough coverage
        char c_base = 'N';
        if(symbol_counts[idx].count >= min_call_coverage)
            c_base = symbol_counts[idx].symbol;
        initial_consensus.append(1, c_base);
    }

    // Copy the final consensus, truncating leading and trailing 'N's
    // If there is an internal no call, the empty string is return
    printf("Initial: %s\n", initial_consensus.c_str());

    std::string consensus;
    bool found_nocall_after_leader = false;
    for(size_t i = 0; i < initial_consensus.size(); ++i)
    {
        if(initial_consensus[i] == 'N')
        {
            // If we have written any bases to the consensus
            // This N represents the trailing sequence or an internal N
            if(!consensus.empty())
                found_nocall_after_leader = true;
        }
        else
        {
            // We have a real base after an internal N
            // return nothing
            if(found_nocall_after_leader)
                return "";
            else
                consensus.push_back(initial_consensus[i]);
        }
    }
    printf("Consens: %s\n", consensus.c_str());
    return consensus;
}
