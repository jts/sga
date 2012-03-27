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
#include "SGVisitors.h"
#include "Profiler.h"
#include "HapgenUtil.h"
#include "multiple_alignment.h"

//
//
//
OverlapHaplotypeBuilder::OverlapHaplotypeBuilder() : m_minIdentity(90.0f), m_minOverlap(41), m_extend_vertices(0)
{
    m_graph = new StringGraph;
}

//
OverlapHaplotypeBuilder::~OverlapHaplotypeBuilder()
{
    m_seedVertices.clear();
    delete m_graph;
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
    
    if(reads.size() > 100)
        return HBRC_OK;

    printf("Starting processing with %zu reads\n", reads.size());
    StringVector ordered_reads;
    addInitialReadsToGraph(reads);
    extendGraph();
    cleanDuplicates();
    extendGraph();
    cleanDuplicates();
    std::string consensus = findHaplotypes();
    if(!consensus.empty())
        out_haplotypes.push_back(consensus);
//    writeGraph("extend-1.dot");
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
}

// Order the list of reads based on the fact that they all share the same kmer
void OverlapHaplotypeBuilder::addInitialReadsToGraph(const StringVector& reads)
{
    // Add reads to graph
    StringVector read_names;
    for(size_t i = 0; i < reads.size(); ++i)
    {
        std::stringstream id_ss;
        id_ss << "base-" << i;
        read_names.push_back(id_ss.str());

        // Add vertex
        printf("OR %s: %s\n", read_names[i].c_str(), reads[i].c_str());
    
        // Save the pointer to the seed vertex, which contains the kmer we started the assembly from
        m_seedVertices.push_back(addVertexToGraph(id_ss.str(), reads[i]));
    }

    // Compute all pairwise overlaps and add edges to the graph
    for(size_t i = 0; i < reads.size(); ++i)
    {
        for(size_t j = i + 1; j < reads.size(); ++j)
        {
            SequenceOverlap overlap = Overlapper::computeOverlap(reads[i], reads[j], ungapped_params);
            bool is_contain_0 = overlap.match[0].start != 0 && overlap.match[0].end != (int)reads[i].size() - 1;
            bool is_contain_1 = overlap.match[1].start != 0 && overlap.match[1].end != (int)reads[j].size() - 1;
            if(overlap.getPercentIdentity() >= m_minIdentity && overlap.getOverlapLength() >= m_minOverlap && !is_contain_0 && !is_contain_1)
            {
                // Make an sga overlap object from the match
                Overlap o(read_names[i], overlap.match[0].start, overlap.match[0].end, reads[i].size(),
                          read_names[j], overlap.match[1].start, overlap.match[1].end, reads[j].size(), false, 0);

                overlap.printAlignment(reads[i], reads[j]);
                SGAlgorithms::createEdgesFromOverlap(m_graph, o, true);
            }
        }
    }

    // Remove duplicate reads and contaiments
    cleanDuplicates();
}

//
Vertex* OverlapHaplotypeBuilder::addVertexToGraph(const std::string& id, const std::string& sequence)
{
    Vertex* pVertex = new(m_graph->getVertexAllocator()) Vertex(id, sequence);
    m_graph->addVertex(pVertex);
    m_used_sequences.insert(sequence);
    return pVertex;
}

// Extend the tips of the graph. Returns a vector
// of pointers to the newly added vertices.
VertexPtrVec OverlapHaplotypeBuilder::extendGraph()
{
    VertexPtrVec added_vertices;
    VertexPtrVec vertices = m_graph->getAllVertices();
    for(size_t i = 0; i < vertices.size(); ++i)
    {
        Vertex* vertex = vertices[i];

        // Check if this vertex needs extension
        size_t as_count = vertex->countEdges(ED_ANTISENSE);
        size_t s_count = vertex->countEdges(ED_SENSE);
        
        // Vertex is bidirectionally connected or not connected at all, skip
        if((as_count > 0 && s_count > 0) || (as_count == 0 && s_count == 0))
            continue;

        EdgeDir dir = as_count == 0 ? ED_ANTISENSE : ED_SENSE;

        printf("Extending %s in dir %d\n", vertex->getID().c_str(), dir);

        // Get overlapping reads using kmers
        StringVector overlapping_reads = getOverlappingReads(vertex->getSeq().toString());
        removeUsedSequences(&overlapping_reads);
        printf("Has %zu overlaps\n", overlapping_reads.size());

        for(size_t j = 0; j < overlapping_reads.size(); ++j)
        {
            Vertex* added = addExtensionVertexToGraph(vertex, overlapping_reads[j]);
            added_vertices.push_back(added);
            bool is_join = isVertexJoinNode(added);

            // Add the vertex to the join list
            // It is added to the opposite direction because we set the direction FROM the node, not TO
            if(is_join)
                m_joinVertices[dir].push_back(added);
            printf("Is %s a join? %d\n", added->getID().c_str(), is_join);
        }
    }

    return added_vertices;
}

//
Vertex* OverlapHaplotypeBuilder::addExtensionVertexToGraph(Vertex* source, const std::string& sequence)
{
    // Add the new vertex to the graph
    std::stringstream id_ss;
    id_ss << "extend-" << m_extend_vertices++;
    Vertex* vertex = addVertexToGraph(id_ss.str(), sequence);

    // Attempt to overlap against all the vertices connected to the source
    // vertex. Potentially slow!
    VertexPtrVec overlap_targets(1, source);
    EdgePtrVec edges = source->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
        overlap_targets.push_back(edges[i]->getEnd());

    for(size_t i = 0; i < overlap_targets.size(); ++i)
    {
        std::string target_sequence = overlap_targets[i]->getSeq().toString();
        SequenceOverlap overlap = Overlapper::computeOverlap(target_sequence, sequence, ungapped_params);
        bool is_contain_0 = overlap.match[0].start != 0 && overlap.match[0].end != (int)target_sequence.size() - 1;
        bool is_contain_1 = overlap.match[1].start != 0 && overlap.match[1].end != (int)sequence.size() - 1;
        if(overlap.getPercentIdentity() >= m_minIdentity && overlap.getOverlapLength() >= m_minOverlap && !is_contain_0 && !is_contain_1)
        {
            // Make an sga overlap object from the match
            Overlap o(overlap_targets[i]->getID(), overlap.match[0].start, overlap.match[0].end, target_sequence.size(),
                      id_ss.str(), overlap.match[1].start, overlap.match[1].end, sequence.size(), false, 0);

            SGAlgorithms::createEdgesFromOverlap(m_graph, o, true);
        }        
    }
    return vertex;
}

//
bool OverlapHaplotypeBuilder::isVertexJoinNode(const Vertex* vertex) const
{
    assert(vertex != NULL);
    std::string sequence = vertex->getSeq().toString();
    size_t nk = sequence.size() - m_parameters.kmer + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer = sequence.substr(i, m_parameters.kmer);
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer, 
                                                                             m_parameters.pBaseBWT, 
                                                                             m_parameters.pBaseBWTCache);

        if(base_count > 0)
            return true;
    }
    return false;
}

//
std::string OverlapHaplotypeBuilder::findHaplotypes()
{

    // Remove transitive edges from the graph
    SGTransitiveReductionVisitor trVisit;
    m_graph->visit(trVisit);

    // Search between each pair of left/right join vertices
    for(size_t i = 0; i < m_joinVertices[ED_ANTISENSE].size(); ++i)
    {
        for(size_t j = 0; j < m_joinVertices[ED_SENSE].size(); ++j)
        {
            Vertex* left_join = m_joinVertices[ED_ANTISENSE].back();//[i];
            Vertex* right_join = m_joinVertices[ED_SENSE].back();

            SGWalkVector walks;
            /*bool success =*/ SGSearch::findWalks(left_join, right_join, ED_SENSE,
                                                   1000, 10000, true, walks);

            printf("Found %zu walks\n", walks.size());
            if(!walks.empty())
            {
                MultipleAlignment multiple_alignment = buildMultipleAlignmentFromWalk(walks.front());
                multiple_alignment.print(200);
                std::string consensus = getConsensus(&multiple_alignment, 5, 1000);
                return consensus;
            }
        }
    }

    return "";
}

//
void OverlapHaplotypeBuilder::removeUsedSequences(StringVector* sequences)
{
    StringVector out;
    for(size_t i = 0; i < sequences->size(); ++i)
    {
        if(m_used_sequences.find(sequences->at(i)) == m_used_sequences.end())
            out.push_back(sequences->at(i));
    }

    sequences->swap(out);
}

//
void OverlapHaplotypeBuilder::cleanDuplicates()
{
    SGDuplicateVisitor duplicateVisit;
    m_graph->visit(duplicateVisit);
    
    // enormous hack
    m_graph->setContainmentFlag(false);
}


//
void OverlapHaplotypeBuilder::writeGraph(const std::string& filename) const
{
    SGTransitiveReductionVisitor trVisit;
    m_graph->visit(trVisit);
    m_graph->writeDot(filename);
}

//
MultipleAlignment OverlapHaplotypeBuilder::buildMultipleAlignmentFromWalk(const SGWalk& walk)
{
    // Insert first sequence
    assert(walk.getNumVertices() != 0);
    MultipleAlignment multiple_alignment;
    multiple_alignment.addBaseSequence(walk.getStartVertex()->getID(), walk.getStartVertex()->getSeq().toString(), "");

    // Insert remaining sequences
    for(size_t i = 1; i < walk.getNumVertices(); ++i)
    {
        std::string s0 = walk.getVertex(i-1)->getSeq().toString();
        std::string s1 = walk.getVertex(i)->getSeq().toString();
        std::string i1 = walk.getVertex(i)->getID();

        SequenceOverlap overlap = Overlapper::computeOverlap(s0, s1, ungapped_params);
        int size_0 = s0.size();
        int size_1 = s1.size();
        bool is_contain_0 = overlap.match[0].start != 0 && overlap.match[0].end != size_0 - 1;
        bool is_contain_1 = overlap.match[1].start != 0 && overlap.match[1].end != size_1 - 1;
        assert(!is_contain_0 && !is_contain_1);
        multiple_alignment.addExtension(i1, s1, "", overlap);
    }

    return multiple_alignment;
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

// Struct to hold a partial match in the FM-index
// The position field is the location in the query sequence of this kmer.
// The index field is an index into the BWT. 
// The is_reverse flag indicates the strand of the partial match
struct KmerMatch
{
    int64_t position:16;
    int64_t index:47;
    int64_t is_reverse:1;

    friend bool operator<(const KmerMatch& a, const KmerMatch& b)
    {
        if(a.index == b.index)
            return a.is_reverse < b.is_reverse;
        else
            return a.index < b.index;
    }

    friend bool operator==(const KmerMatch& a, const KmerMatch& b)
    {
        return a.index == b.index && a.is_reverse == b.is_reverse;
    }
};

// Return a hash key for a KmerMatch
struct KmerMatchKey
{
    size_t operator()(const KmerMatch& a) const { return a.index; }
};

typedef std::set<KmerMatch> KmerMatchSet;
typedef HashMap<KmerMatch, bool, KmerMatchKey> KmerMatchMap;

StringVector OverlapHaplotypeBuilder::getOverlappingReads(const std::string& sequence) const
{
    WARN_ONCE("TODO: Refactor the kmer-based overlap code in OverlapHaplotyBuilder::getOverlappingReads")
    StringVector reads;
    int64_t max_interval_size = 500;

    // Use the FM-index to look up intervals for each kmer of the read. Each index
    // in the interval is stored individually in the KmerMatchMap. We then
    // backtrack to map these kmer indices to read IDs. As reads can share
    // multiple kmers, we use the map to avoid redundant lookups.
    // There is likely a faster algorithm which performs direct decompression
    // of the read sequences without having to expand the intervals to individual
    // indices. The current algorithm suffices for now.
    KmerMatchMap prematchMap;
    size_t num_kmers = sequence.size() - m_parameters.kmer + 1;
    for(size_t i = 0; i < num_kmers; ++i)
    {
        std::string kmer = sequence.substr(i, m_parameters.kmer);
        BWTInterval interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, kmer);
        if(interval.isValid() && interval.size() < max_interval_size) 
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
            {
                KmerMatch match = { i, j, false };
                prematchMap.insert(std::make_pair(match, false));
            }
        }

        kmer = reverseComplement(kmer);
        interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, kmer);
        if(interval.isValid() && interval.size() < max_interval_size) 
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
            {
                KmerMatch match = { i, j, true };
                prematchMap.insert(std::make_pair(match, false));
            }
        }
    }

    // Backtrack through the kmer indices to turn them into read indices.
    // This mirrors the calcSA function in SampledSuffixArray except we mark each entry
    // as visited once it is processed.
    KmerMatchSet matches;
    for(KmerMatchMap::iterator iter = prematchMap.begin(); iter != prematchMap.end(); ++iter)
    {
        // This index has been visited
        if(iter->second)
            continue;

        // Mark this as visited
        iter->second = true;

        // Backtrack the index until we hit the starting symbol
        KmerMatch out_match = iter->first;
        while(1) 
        {
            char b = m_parameters.pVariantBWT->getChar(out_match.index);
            out_match.index = m_parameters.pVariantBWT->getPC(b) + m_parameters.pVariantBWT->getOcc(b, out_match.index - 1);

            // Check if the hash indicates we have visited this index. If so, stop the backtrack
            KmerMatchMap::iterator find_iter = prematchMap.find(out_match);
            if(find_iter != prematchMap.end())
            {
                // We have processed this index already
                if(find_iter->second)
                    break;
                else
                    find_iter->second = true;
            }

            if(b == '$')
            {
                // We've found the lexicographic index for this read. Turn it into a proper ID
                out_match.index = m_parameters.pVariantSSA->lookupLexoRank(out_match.index);
                matches.insert(out_match);
                break;
            }
        }
    }

    // Perform the actual extraction from the index
    for(KmerMatchSet::iterator iter = matches.begin(); iter != matches.end(); ++iter)
    {
        std::string match_sequence = BWTAlgorithms::extractString(m_parameters.pVariantBWT, iter->index);
        if(iter->is_reverse)
            match_sequence = reverseComplement(match_sequence);
        reads.push_back(match_sequence);
    }

    return reads;
}

