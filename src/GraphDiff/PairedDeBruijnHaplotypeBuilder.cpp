///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// PairedDeBruijnHaplotypeBuilder - Construct candidate
// haplotypes from a pair of k-mer seeds.
//
#include "PairedDeBruijnHaplotypeBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Profiler.h"
#include "HapgenUtil.h"
#include "multiple_alignment.h"
#include "KmerOverlaps.h"
#include "Verbosity.h"

//
//
//
PairedDeBruijnHaplotypeBuilder::PairedDeBruijnHaplotypeBuilder(const GraphCompareParameters& params) : m_parameters(params)
{

}

//
PairedDeBruijnHaplotypeBuilder::~PairedDeBruijnHaplotypeBuilder()
{

}

//
void PairedDeBruijnHaplotypeBuilder::setInitialHaplotype(const std::string& str)
{
    m_startingKmer = str;
}

// Run the bubble construction process
HaplotypeBuilderReturnCode PairedDeBruijnHaplotypeBuilder::run(StringVector& out_haplotypes)
{
    PROFILE_FUNC("GraphCompare::buildVariantStringGraph")
    assert(!m_startingKmer.empty());
    
    std::map<std::string, int> kmerCountMap;

    // We search until we find the first common vertex in each direction
    //size_t MIN_TARGET_COUNT = m_parameters.bReferenceMode ? 1 : 2;
    size_t MAX_ITERATIONS = 50000;
    size_t MAX_SIMULTANEOUS_BRANCHES = 100;
    size_t MAX_TOTAL_BRANCHES = 500;

    // Tracking stats
    size_t max_simul_branches_used = 0;
    size_t total_branches = 0;
    size_t iterations = 0;

    // Select kmers from the pairs of reads containing the variant kmers
    // We search the graph until we find these targets
    std::set<std::string> target_set;
    selectTargetKmers(m_startingKmer, true, target_set);
    selectTargetKmers(reverseComplement(m_startingKmer), false, target_set);
    
    if(Verbosity::Instance().getPrintLevel() > 3)
        printf("PairedDBGHaplotype: found %zu targets\n", target_set.size());

    if(target_set.empty())
        return HBRC_OK;

    // Initialize the graph
    StringGraph* pGraph = new StringGraph;
    BuilderExtensionQueue queue;

    Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(m_startingKmer, m_startingKmer);
    pVertex->setColor(GC_BLACK);
    pGraph->addVertex(pVertex);

    // Add the vertex to the extension queue
    queue.push(BuilderExtensionNode(pVertex, ED_SENSE));
    queue.push(BuilderExtensionNode(pVertex, ED_ANTISENSE));

    std::vector<Vertex*> sense_join_vector;
    std::vector<Vertex*> antisense_join_vector;

    bool extension_allowed[ED_COUNT] = { true, true };

    // Perform the extension. The while conditions are heuristics to avoid searching
    // the graph too much 
    while(!queue.empty() && iterations++ < MAX_ITERATIONS && queue.size() < MAX_SIMULTANEOUS_BRANCHES && total_branches < MAX_TOTAL_BRANCHES)
    {
        if(queue.size() > max_simul_branches_used)
            max_simul_branches_used = queue.size();

        BuilderExtensionNode curr = queue.front();
        queue.pop();

        if(!extension_allowed[curr.direction])
            continue;

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_parameters.variantIndex.pBWT, curr.direction);

        std::string extensionsUsed;
        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            bool acceptExt = count >= m_parameters.minDBGCount;
            if(!acceptExt)
                continue;

            extensionsUsed.push_back(b);
            std::string newStr = VariationBuilderCommon::makeDeBruijnVertex(vertStr, b, curr.direction);
            kmerCountMap[newStr] = count;

            // Create the new vertex and edge in the graph
            // Skip if the vertex already exists
            if(pGraph->getVertex(newStr) != NULL)
                continue;
            
            // Allocate the new vertex and add it to the graph
            Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(newStr, newStr);
            pVertex->setColor(GC_BLACK);
            pGraph->addVertex(pVertex);

            // Add edges
            VariationBuilderCommon::addSameStrandDeBruijnEdges(pGraph, curr.pVertex, pVertex, curr.direction);
            
            if(target_set.find(newStr) != target_set.end())
            {
                if(curr.direction == ED_SENSE)
                    sense_join_vector.push_back(pVertex);
                else
                    antisense_join_vector.push_back(pVertex);
                extension_allowed[curr.direction] = false;
            }
            else
            {
                // Add the vertex to the extension queue
                queue.push(BuilderExtensionNode(pVertex, curr.direction));
            }
        }
        
        // Update the total number of times we branches the search
        if(!extensionsUsed.empty())
            total_branches += extensionsUsed.size() - 1;
    }
    
    if(Verbosity::Instance().getPrintLevel() > 2)
    {
        printf("[PairedDBG] iterations: %zu total_branches: %zu left found? %d right found? %d\n", 
                 iterations, total_branches, !extension_allowed[ED_ANTISENSE], !extension_allowed[ED_SENSE]);
    }

    // If the graph construction was successful, walk the graph
    // between the endpoints to make a string
    // Generate haplotypes between every pair of antisense/sense join vertices
    for(size_t i = 0; i < antisense_join_vector.size(); ++i) {
        for(size_t j = 0; j < sense_join_vector.size(); ++j) {
            SGWalkVector outWalks;
            SGSearch::findWalks(antisense_join_vector[i],
                                sense_join_vector[j],
                                ED_SENSE,
                                100000, // max distance to search
                                10000, // max nodes to search
                                true, // exhaustive search
                                outWalks);

            for(size_t k = 0; k < outWalks.size(); ++k)
                out_haplotypes.push_back(outWalks[k].getString(SGWT_START_TO_END));
        }
    }

    if(!out_haplotypes.empty())
    {
        std::cout << "haplotypes found!\n";
        for(size_t i = 0; i < out_haplotypes.size(); i++)
        {
            printf(">HID:%zu\n%s\n", i, out_haplotypes[i].c_str());
        }
    }

    delete pGraph;
    return HBRC_OK;
}
 
//
void PairedDeBruijnHaplotypeBuilder::selectTargetKmers(const std::string& kmer, bool rc_targets, std::set<std::string>& target_set) const
{
    size_t k = kmer.size();
    std::vector<size_t> ids = getReadIDs(kmer);

    // Iterate over the pairs of the reads containing the input kmer
    // and return the first reference kmer in any pair
    for(size_t i = 0; i < ids.size(); i++)
    {
        size_t pair_id = ids[i] % 2 == 0 ? ids[i] + 1 : ids[i] - 1;
        std::string pair = BWTAlgorithms::extractString(m_parameters.variantIndex.pBWT, pair_id);
        printf(">PID:%zu\n%s\n", pair_id, pair.c_str());

        for(size_t j = 0; j < pair.size() - k + 1; ++j)
        {
            std::string test_kmer = pair.substr(j, k);
            size_t ref_count = BWTAlgorithms::countSequenceOccurrences(test_kmer, m_parameters.referenceIndex);
            size_t base_count = BWTAlgorithms::countSequenceOccurrences(test_kmer, m_parameters.baseIndex);
            if(ref_count == 1 && base_count > 0)
                target_set.insert(rc_targets ? reverseComplement(test_kmer) : test_kmer);
        }
    }
}

//       
std::vector<size_t> PairedDeBruijnHaplotypeBuilder::getReadIDs(const std::string& kmer) const
{
    BWTInterval interval = BWTAlgorithms::findInterval(m_parameters.variantIndex, kmer);
    std::cout << "Kmer: " << kmer << "\n";
    std::vector<size_t> out;
    for(int64_t i = interval.lower; i <= interval.upper; ++i)
    {
        out.push_back(m_parameters.variantIndex.pSSA->calcSA(i, m_parameters.variantIndex.pBWT).getID());
        std::cout << "ID: " << out.back() << "\n";
    }
    return out;
}
