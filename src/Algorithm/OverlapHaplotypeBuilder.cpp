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
#include "KmerOverlaps.h"

//
//
//
OverlapHaplotypeBuilder::OverlapHaplotypeBuilder(const GraphCompareParameters& params) : m_parameters(params), m_numReads(0)
{
    assert(m_parameters.minOverlap > 25);

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
    if(!buildInitialGraph(reads))
        return HBRC_OK;

    bool done = false;
    int MAX_ROUNDS = 4;
    int round = 0;
    while(!done) 
    {
        // Extend the graph by finding new overlaps for the reads at the tips
        extendGraph();

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
        
        size_t num_vertices = m_graph->getNumVertices();
        size_t num_tips = findTips().size();

        done = round++ >= MAX_ROUNDS || !out_haplotypes.empty() || num_tips > 5;

        printf("graph size: %zu tips: %zu\n", num_vertices, num_tips);
    }

    m_graph->writeDot("final.dot");

    return HBRC_OK;
}

//
bool OverlapHaplotypeBuilder::buildInitialGraph(const StringVector& reads)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::buildInitialGraph")
    // Compute initial ordering of reads based on the position of the
    // starting kmer sequence. If the starting kmer was corrected out
    // of a read, it is discarded.
    StringVector ordered_reads;
    orderReadsInitial(m_initial_kmer_string, reads, &ordered_reads);

    if(ordered_reads.size() < m_parameters.kmerThreshold)
        return false;

    // DEBUG print MA
    //MultipleAlignment ma = buildMultipleAlignment(ordered_reads);
    //ma.print(200);

    // Insert initial reads into graph
    for(size_t i = 0; i < ordered_reads.size(); ++i)
        insertVertexIntoGraph("seed-", ordered_reads[i]);
    return true;
}

//
void OverlapHaplotypeBuilder::extendGraph()
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::extendGraph")

    // Find tip vertices
    ExtendableTipVector tips = findTips();
    for(size_t i = 0; i < tips.size(); ++i)
    {
        Vertex* x = tips[i].vertex;
        EdgeDir dir = tips[i].direction;

        printf("Vertex %s can be extended\n", x->getID().c_str());

        // Find corrected reads that perfectly overlap this vertex
        StringVector overlapping_reads = getCorrectedOverlaps(x->getSeq().toString());
        printf("Found %zu overlaps\n", overlapping_reads.size());

        // Insert the new reads into the graph
        for(size_t i = 0; i < overlapping_reads.size(); ++i)
        {
            // Determine whether the incoming read is possible join point
            bool is_join = isJoinSequence(overlapping_reads[i]);

            std::stringstream label;
            label << (is_join ? "join-" : "extend-");
            label << (dir == ED_SENSE ? "sense-" : "antisense-");
            insertVertexIntoGraph(label.str(), overlapping_reads[i]);
        }
    }
}

// 
void OverlapHaplotypeBuilder::insertVertexIntoGraph(const std::string& prefix, const std::string& sequence)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::insertVertexIntoGraph")
    // Check if a vertex with this sequence already exists, if so
    // we do not create a new vertex
    if(m_used_reads.find(sequence) != m_used_reads.end())
        return;

    // Create the vertex
    std::stringstream id_ss;
    id_ss << prefix << m_numReads++;
    Vertex* pVertex = new(m_graph->getVertexAllocator()) Vertex(id_ss.str(), sequence);
    m_graph->addVertex(pVertex);

    // Get a list of vertices that have a k-mer match to this sequence.
    SharedVertexKmerVector candidate_vertices = getCandidateOverlaps(pVertex);

    // Align the sequence against candidate reads
    for(size_t i = 0; i < candidate_vertices.size(); ++i) 
    {
        Vertex* existing_vertex = candidate_vertices[i].vertex;

        // Skip self-edge
        if(existing_vertex == pVertex)
            continue;

        std::string existing_sequence = existing_vertex->getSeq().toString();

        // Try to compute the overlap using the shared kmer as a seed match
        SequenceOverlap overlap;
        std::string shared_kmer = sequence.substr(candidate_vertices[i].kmer_index, m_vertex_map_kmer);
        size_t pos_0 = existing_sequence.find(shared_kmer);
        size_t pos_1 = sequence.find(shared_kmer);
        assert(pos_0 != std::string::npos && pos_1 != std::string::npos);

        // Check for secondary occurrences
        if(existing_sequence.find(shared_kmer, pos_0 + 1) != std::string::npos || 
                    sequence.find(shared_kmer, pos_1 + 1) != std::string::npos) 
        {
            // One of the reads has a second occurrence of the kmer. Use
            // the slow overlapper.
            overlap = Overlapper::computeOverlap(existing_sequence, sequence);
        } 
        else 
        {
            // Seed the match using kmer position
            overlap = Overlapper::extendMatch(existing_sequence, sequence, pos_0, pos_1, 1);
        }

        if(overlap.edit_distance == 0 && overlap.getOverlapLength() >= m_parameters.minOverlap)
        {
            // Add an overlap to the graph
            // Translate the sequence overlap struture into an SGA overlap
            Overlap sga_overlap(existing_vertex->getID(), overlap.match[0].start, overlap.match[0].end, existing_sequence.size(),
                                             id_ss.str(), overlap.match[1].start, overlap.match[1].end, sequence.size(), false, 0);
            SGAlgorithms::createEdgesFromOverlap(m_graph, sga_overlap, true);
        }
    }

    // Insert the sequence into the used reads set
    m_used_reads.insert(sequence);

    // Update kmer map too
    updateKmerVertexMap(pVertex);
}

//
SharedVertexKmerVector OverlapHaplotypeBuilder::getCandidateOverlaps(const Vertex* incoming_vertex)
{
    // Get a list of vertices that have a k-mer match to this sequence.
    SharedVertexKmerVector candidate_vertices;
    const std::string sequence = incoming_vertex->getSeq().toString();
    size_t nk = sequence.size() - m_vertex_map_kmer + 1;
    for(size_t i = 0; i < nk; ++i) 
    {
        std::string query_kmer = sequence.substr(i, m_vertex_map_kmer);
        HashMap<std::string, StringVector>::iterator find_iter = m_kmer_vertex_id_cache.find(query_kmer);
        if(find_iter != m_kmer_vertex_id_cache.end())
        {
            // Get the vertices for these ids
            for(StringVector::iterator iter = find_iter->second.begin(); iter != find_iter->second.end(); ++iter)
            {
                Vertex* candidate = m_graph->getVertex(*iter);
                if(candidate != NULL)
                {
                    SharedVertexKmer svk = { i, candidate };
                    candidate_vertices.push_back(svk);
                }
            }
        }
    }

    // Sort the vector by the vertex pointer and remove duplicates
    // Remove duplicates from the candidate list
    std::sort(candidate_vertices.begin(), candidate_vertices.end(), SharedVertexKmer::sortByVertex);
    SharedVertexKmerVector::iterator new_end = std::unique(candidate_vertices.begin(), candidate_vertices.end(), SharedVertexKmer::equalByVertex);
    candidate_vertices.resize(new_end - candidate_vertices.begin());
    return candidate_vertices;
}

// 
void OverlapHaplotypeBuilder::updateKmerVertexMap(const Vertex* incoming_vertex)
{
    const std::string sequence = incoming_vertex->getSeq().toString();
    size_t nk = sequence.size() - m_vertex_map_kmer + 1;

    // Insert kmers into the kmer pointer cache
    for(size_t i = 0; i < nk; ++i) 
    {
        std::string query_kmer = sequence.substr(i, m_vertex_map_kmer);
        m_kmer_vertex_id_cache[query_kmer].push_back(incoming_vertex->getID());
    }
}

// 
void OverlapHaplotypeBuilder::checkWalks(StringVector* walk_strings)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::checkWalks")
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

    printf("%zu walk candidates\n", join_vertices.size());

    //
    for(size_t i = 0; i < join_vertices.size(); ++i)
    {
        for(size_t j = i + 1; j < join_vertices.size(); ++j)
        {
            // Try to find a walk between this pair of join vertices
            SGWalkVector walks;
            SGSearch::findWalks(join_vertices[i], join_vertices[j], ED_SENSE, 2000, 10000, true, walks);

            if(!walks.empty())
            {
                // Check if the walk set is self-contained
                bool self_contained = areWalksValid(walks);
                if(!self_contained)
                    continue;

                // Check how many seed vertices this walk covers
                size_t seeds_covered = countCoveredVertices(seed_vertices, walks);
                //printf("    walk covers %zu of %zu seeds\n", seeds_covered, seed_vertices.size());
                double fraction_covered = (double)seeds_covered / seed_vertices.size();
                if(fraction_covered > 0.5)
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
bool OverlapHaplotypeBuilder::areWalksValid(const SGWalkVector& walks)
{
    // Build ID set
    std::set<Vertex*> vertex_set;
    for(size_t i = 0; i < walks.size(); ++i)
    {
        VertexPtrVec verts = walks[i].getVertices();
        for(size_t j = 0; j < verts.size(); ++j)
            vertex_set.insert(verts[j]);
    }

    // For every internal vertex of each walk, check that it only has edges to vertices in the set
    for(size_t i = 0; i < walks.size(); ++i)
    {
        VertexPtrVec verts = walks[i].getVertices();
        for(size_t j = 1; j < verts.size() - 1; ++j)
        {
            EdgePtrVec edges = verts[j]->getEdges();
            for(size_t k = 0; k < edges.size(); ++k)
            {
                if(vertex_set.find(edges[k]->getEnd()) == vertex_set.end())
                    return false; // edge to a non-contained vertex
            }
        }
    }
    return true;
}

//
ExtendableTipVector OverlapHaplotypeBuilder::findTips() const
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::findTips")
    ExtendableTipVector tips;
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

        // Skip completely disconnected vertices
        if(edge_count[ED_SENSE] == 0 && edge_count[ED_ANTISENSE] == 0)
            continue;

        if(edge_count[ED_SENSE] == 0)
        {
            ExtendableTip tip = { x, ED_SENSE };
            tips.push_back(tip);
        }

        if(edge_count[ED_ANTISENSE] == 0)
        {
            ExtendableTip tip = { x, ED_ANTISENSE };
            tips.push_back(tip);
        }
    }
    return tips;
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
    PROFILE_FUNC("OverlapHaplotypeBuilder::getCorrectedOverlaps")
    // Extract reads that share a short kmer with the input sequence
    size_t k = 41;

    SequenceOverlapPairVector overlap_vector = KmerOverlaps::retrieveMatches(sequence,
                                                                             k,
                                                                             m_parameters.minOverlap,
                                                                             0.95,
                                                                             2,
                                                                             m_parameters.pVariantBWT,
                                                                             m_parameters.pVariantBWTCache,
                                                                             m_parameters.pVariantSSA);

    // Extract the read sequences from the multiple alignment
    StringVector reads;
    for(size_t i = 0; i < overlap_vector.size(); ++i)
        reads.push_back(overlap_vector[i].sequence[1]);
    printf("Extracted %zu reads\n", overlap_vector.size());

    // Correct the reads
    correctReads(&reads);
    assert(reads.size() == overlap_vector.size());
    // Overlap each corrected read with the input sequence
    // Discard any read that does perfectly match with a long enough overlap
    StringVector out_reads;
    for(size_t i = 0; i < reads.size(); ++i) 
    {
        // Skip uncorrected reads
        if(reads[i].empty())
            continue;

        // Recompute the overlap, using the previous overlap as a guide
        const SequenceOverlap& initial_overlap = overlap_vector[i].overlap;
        SequenceOverlap overlap = Overlapper::extendMatch(sequence, reads[i], initial_overlap.match[0].start, initial_overlap.match[1].start, 1);
        if(overlap.edit_distance == 0 && overlap.getOverlapLength() >= m_parameters.minOverlap)
            out_reads.push_back(reads[i]);
    }

    return out_reads;
}

//
void OverlapHaplotypeBuilder::getReadsForKmers(const StringVector& kmer_vector, size_t k, StringVector* reads)
{
    PROFILE_FUNC("OverlapHaplotypeBuilder::getReadsForKmers")
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
    PROFILE_FUNC("OverlapHaplotypeBuilder::correctReads")
    StringVector out_reads;
    for(size_t i = 0; i < reads->size(); ++i)
    {
        // Check if this read has been corrected before
        std::string corrected_sequence;
        HashMap<std::string, std::string>::iterator find_iter = m_correction_cache.find(reads->at(i));
        if(find_iter != m_correction_cache.end())
        {
            corrected_sequence = find_iter->second;
        }
        else
        {
            // Perform correction
            DNAString dna(reads->at(i));
            SeqRecord record = { "null", dna, "" };
            SequenceWorkItem wi(0, record);
            ErrorCorrectResult r = m_corrector->correct(wi);
            if(r.kmerQC || r.overlapQC)
                corrected_sequence = r.correctSequence.toString();
            
            // Insert the sequence into the cache
            // We allow it to be empty 
            m_correction_cache.insert(std::make_pair(reads->at(i), corrected_sequence));
        }

        // We allow the sequence to be null if the read could not be corrected
        out_reads.push_back(corrected_sequence);
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

    // 
    if(read_kmer_vector.empty())
        return;

    // Sort the vector by the second element
    std::sort(read_kmer_vector.begin(), read_kmer_vector.end(), sortPairSecondAscending);

    // Insert the reads into the ordered list
    for(size_t i = 0; i < read_kmer_vector.size(); ++i)
        ordered_vector->push_back(reads[read_kmer_vector[i].first]);
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

