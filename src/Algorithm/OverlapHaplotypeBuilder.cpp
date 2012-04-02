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
OverlapHaplotypeBuilder::OverlapHaplotypeBuilder(const GraphCompareParameters& params) : m_parameters(params), m_numReads(0)
{
    ErrorCorrectParameters correction_params;
    correction_params.pOverlapper = NULL;
    correction_params.pBWT = m_parameters.pVariantBWT;
    correction_params.pSSA = m_parameters.pVariantSSA;
    assert(m_parameters.pVariantBWTCache != NULL);
    correction_params.pIntervalCache = m_parameters.pVariantBWTCache;
    correction_params.algorithm = ECA_KMER;

    // Overlap-based corrector params
    correction_params.minOverlap = 51;
    correction_params.numOverlapRounds = 1;
    correction_params.minIdentity = 95.0f;
    correction_params.conflictCutoff = 5;
    correction_params.depthFilter = 100;

    // k-mer based corrector params
    correction_params.numKmerRounds = 10;
    correction_params.kmerLength = 51;

    // output options
    correction_params.printOverlaps = true;

    m_graph = new StringGraph;
    m_corrector = new ErrorCorrectProcess(correction_params);
}

//
OverlapHaplotypeBuilder::~OverlapHaplotypeBuilder()
{
    delete m_corrector;
    delete m_graph;
}

// The source string is the string the bubble starts from
void OverlapHaplotypeBuilder::setInitialHaplotype(const std::string& sequence)
{
    printf("\n\n***** Starting new haplotype %s\n", sequence.c_str());
    m_initial_kmer_string = sequence;
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
    getReadsForKmers(query_kmers, m_parameters.kmer, &reads);
    printf("Found %zu initial reads\n", reads.size());
    // Correct sequencing errors in the reads
    correctReads(&reads);
    printf("Corrected %zu initial reads\n", reads.size());
    
    // Start the graph using the corrected reads that contain the initial kmer
    buildInitialGraph(reads);

    // Extend the graph by finding new overlaps for the reads at the tips
    extendGraph();
    extendGraph();
    extendGraph();

    m_graph->writeASQG("hapgraph.asqg");
    
    // Clean up the graph
    SGIdenticalRemoveVisitor dupVisit;
    m_graph->visit(dupVisit);

    // hack
    WARN_ONCE("hacked containment flag");
    m_graph->setContainmentFlag(false);

    SGTransitiveReductionVisitor trVisit;
    m_graph->visit(trVisit);

    // Check for walks that cover all seed vertices
    checkWalks(&out_haplotypes);
    m_graph->writeDot("hapgraph.dot");
    return HBRC_OK;
}

//
void OverlapHaplotypeBuilder::buildInitialGraph(const StringVector& reads)
{
    // Compute initial ordering of reads based on the position of the
    // starting kmer sequence. If the starting kmer was corrected out
    // of a read, it is discarded.
    StringVector ordered_reads;
    orderReadsInitial(m_initial_kmer_string, reads, &ordered_reads);

    if(ordered_reads.empty())
        return;

    // DEBUG print MA
    MultipleAlignment ma = buildMultipleAlignment(ordered_reads);
    ma.print(200);

    // Insert the first read into the graph
    std::stringstream id_ss;
    id_ss << "seed-" << m_numReads++;
    Vertex* pVertex = new(m_graph->getVertexAllocator()) Vertex(id_ss.str(), ordered_reads.front());
    m_graph->addVertex(pVertex);

    // Insert all other reads into the graph
    for(size_t i = 1; i < ordered_reads.size(); ++i)
        insertVertexIntoGraph("seed-", ordered_reads[i]);

}

//
void OverlapHaplotypeBuilder::extendGraph()
{

    // Find tip vertices
    VertexPtrVec tips = findTips();
    for(size_t i = 0; i < tips.size(); ++i)
    {
        printf("Vertex %s can be extended\n", tips[i]->getID().c_str());

        // Find corrected reads that perfectly overlap this vertex
        StringVector overlapping_reads = getCorrectedOverlaps(tips[i]->getSeq().toString());
        printf("Found %zu overlaps\n", overlapping_reads.size());

        // Insert the new reads into the graph
        for(size_t i = 0; i < overlapping_reads.size(); ++i)
        {
            // Determine whether the incoming read is possible join point
            bool is_join = isJoinSequence(overlapping_reads[i]);
            insertVertexIntoGraph(is_join ? "join-" : "extend-", overlapping_reads[i]);
        }
    }
}

// 
void OverlapHaplotypeBuilder::insertVertexIntoGraph(const std::string& prefix, const std::string& sequence)
{
    // Check if a vertex with this sequence already exists, if so
    // we do not create a new vertex
    if(m_used_reads.find(sequence) != m_used_reads.end())
        return;

    // Create the vertex
    std::stringstream id_ss;
    id_ss << prefix << m_numReads++;
    Vertex* pVertex = new(m_graph->getVertexAllocator()) Vertex(id_ss.str(), sequence);
    m_graph->addVertex(pVertex);

    // Slow function which tests the incoming read against everything in the graph
    VertexPtrVec vertices = m_graph->getAllVertices();
    for(size_t i = 0; i < vertices.size(); ++i) 
    {
        // Skip self-edge
        if(vertices[i] == pVertex)
            continue;

        // TODO: replace this with fast function
        std::string v_sequence = vertices[i]->getSeq().toString();
        SequenceOverlap overlap = Overlapper::computeOverlap(v_sequence, sequence, ungapped_params);
        if(overlap.edit_distance == 0 && overlap.getOverlapLength() >= MIN_OVERLAP)
        {
            // Add an overlap to the graph
            // Translate the sequence overlap struture into an SGA overlap
            Overlap sga_overlap(vertices[i]->getID(), overlap.match[0].start, overlap.match[0].end, v_sequence.size(),
                                         id_ss.str(), overlap.match[1].start, overlap.match[1].end, sequence.size(), false, 0);
            SGAlgorithms::createEdgesFromOverlap(m_graph, sga_overlap, true);
        }
    }

    // Insert the sequence into the used reads set
    m_used_reads.insert(sequence);
}

// 
void OverlapHaplotypeBuilder::checkWalks(StringVector* walk_strings)
{
    // Make a vector of the join sequences
    VertexPtrVec seed_vertices;
    VertexPtrVec join_vertices;
    VertexPtrVec vertices = m_graph->getAllVertices();
    for(size_t i = 0; i < vertices.size(); ++i)
    {
        if(vertices[i]->getID().find("join") != std::string::npos)
            join_vertices.push_back(vertices[i]);
        if(vertices[i]->getID().find("seed") != std::string::npos)
            seed_vertices.push_back(vertices[i]);
    }

    printf("Found %zu potential joins\n", join_vertices.size());
    for(size_t i = 0; i < join_vertices.size(); ++i)
    {
        for(size_t j = i + 1; j < join_vertices.size(); ++j)
        {
            // Try to find a walk between this pair of join vertices
            SGWalkVector walks;
            SGSearch::findWalks(join_vertices[i], join_vertices[j], ED_SENSE, 2000, 10000, true, walks);

            if(!walks.empty())
            {
                printf("Joined walks found\n");

                // Check how many seed vertices this walk covers
                size_t seeds_covered = countCoveredVertices(seed_vertices, walks);
                printf("    walk covers %zu of %zu seeds\n", seeds_covered, seed_vertices.size());
                if(seeds_covered == seed_vertices.size())
                {
                    // Add the sequence of the completed walks to the output vector
                    for(size_t i = 0; i < walks.size(); ++i)
                        walk_strings->push_back(walks[i].getString(SGWT_START_TO_END));
                    return;
                }
            }
        }
    }
}

//
VertexPtrVec OverlapHaplotypeBuilder::findTips() const
{
    VertexPtrVec out_vertices;

    VertexPtrVec vertices = m_graph->getAllVertices();
    for(size_t i = 0; i < vertices.size(); ++i)
    {
        size_t edge_count[ED_COUNT];
        edge_count[ED_SENSE] = 0;
        edge_count[ED_ANTISENSE] = 0;

        Vertex* x = vertices[i];
        EdgePtrVec x_edges = x->getEdges();
        for(size_t j = 0; j < x_edges.size(); ++j)
        {
            Edge* xy = x_edges[j];

            // Do not count edges to duplicate vertices
            Overlap o = xy->getOverlap();
            if(!o.isContainment())
                edge_count[xy->getDir()]++;
        }

        if(edge_count[ED_SENSE] == 0 || edge_count[ED_ANTISENSE] == 0)
            out_vertices.push_back(x);
    }
    return out_vertices;
}

bool OverlapHaplotypeBuilder::isJoinSequence(const std::string& sequence)
{
    // Check if this sequence forms a complete k-mer path in the base/reference sequence
    size_t k = m_parameters.kmer;
    if(sequence.size() < k)
        return false;

    size_t nk = sequence.size() - k + 1;
    for(size_t i = 0; i < nk; ++i)
    {
        std::string kmer_seq = sequence.substr(i, k);
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer_seq, 
                                                                             m_parameters.pBaseBWT, 
                                                                             m_parameters.pBaseBWTCache);

        if(base_count == 0)
            return false;
    }

    return true;
}

//
StringVector OverlapHaplotypeBuilder::getCorrectedOverlaps(const std::string& sequence)
{
    // Extract reads that share a short kmer with the input sequence
    size_t k = 31;
    size_t nk = sequence.size() - k + 1;
    StringVector query_kmers;
    for(size_t i = 0; i < nk; ++i)
        query_kmers.push_back(sequence.substr(i, k));

    StringVector reads;
    getReadsForKmers(query_kmers, k, &reads);
    printf("    %zu initial reads\n", reads.size());
    // Correct the reads
    correctReads(&reads);

    // Overlap each corrected read with the input sequence
    // Discard any read that does perfectly match with a long enough overlap
    StringVector out_reads;
    for(size_t i = 0; i < reads.size(); ++i) 
    {
        SequenceOverlap overlap = Overlapper::computeOverlap(sequence, reads[i], ungapped_params);
        if(overlap.edit_distance == 0 && overlap.getOverlapLength() >= MIN_OVERLAP)
            out_reads.push_back(reads[i]);
    }

    return out_reads;
}

//
void OverlapHaplotypeBuilder::getReadsForKmers(const StringVector& kmer_vector, size_t k, StringVector* reads)
{
    SeqItemVector fwd_si;
    SeqItemVector rev_si;

    // Forward reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, k, false, 100000, &fwd_si, NULL);

    // Reverse reads
    HapgenUtil::extractHaplotypeReads(kmer_vector, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, 
                                      m_parameters.pVariantSSA, k, true, 100000, &rev_si, NULL);

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

//
size_t OverlapHaplotypeBuilder::countCoveredVertices(const VertexPtrVec& vertices, const SGWalkVector& walks)
{
    size_t count = 0;
    std::set<Vertex*> walk_vertices;
    for(size_t i = 0; i < walks.size(); ++i)
    {
        for(size_t j = 0; j < walks[i].getNumVertices(); ++j)
            walk_vertices.insert(walks[i].getVertex(j));
    }

    for(size_t i = 0; i < vertices.size(); ++i)
    {
        if(walk_vertices.find(vertices[i]) != walk_vertices.end())
            count += 1;
    }

    return count;
}

//
void OverlapHaplotypeBuilder::correctReads(StringVector* reads)
{
    StringVector out_reads;
    for(size_t i = 0; i < reads->size(); ++i)
    {
        DNAString dna(reads->at(i));
        SeqRecord record = { "null", dna, "" };
        SequenceWorkItem wi(0, record);
        ErrorCorrectResult r = m_corrector->correct(wi);
        if(r.kmerQC || r.overlapQC)
            out_reads.push_back(r.correctSequence.toString());
    }
    reads->swap(out_reads);
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
        if(pos != std::string::npos)
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
        multiple_alignment.addExtension("seq", *read_iterator, "", overlap);
        prev_iterator = read_iterator;
        read_iterator++;
    }

    return multiple_alignment;
}

