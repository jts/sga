///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// DeBruijnHaplotypeBuilder - Construct candidate
// haplotypes from a pair of k-mer seeds.
//
#include "DeBruijnHaplotypeBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "SGVisitors.h"
#include "Profiler.h"
#include "HapgenUtil.h"
#include "multiple_alignment.h"
#include "KmerOverlaps.h"
#include "CorrectionThresholds.h"

//
//
//
DeBruijnHaplotypeBuilder::DeBruijnHaplotypeBuilder(const GraphCompareParameters& params) : m_parameters(params)
{

}

//
DeBruijnHaplotypeBuilder::~DeBruijnHaplotypeBuilder()
{

}

//
void DeBruijnHaplotypeBuilder::setInitialHaplotype(const std::string& str)
{
    m_startingKmer = str;
}

// Run the bubble construction process
HaplotypeBuilderReturnCode DeBruijnHaplotypeBuilder::run(StringVector& out_haplotypes)
{
    PROFILE_FUNC("GraphCompare::buildVariantStringGraph")
    assert(!m_startingKmer.empty());

    std::map<std::string, int> kmerCountMap;

    // We search until we find the first common vertex in each direction
    size_t MIN_TARGET_COUNT = m_parameters.bReferenceMode ? 1 : 2;
    size_t MAX_ITERATIONS = 2000;
    size_t MAX_SIMULTANEOUS_BRANCHES = 40;
    size_t MAX_TOTAL_BRANCHES = 50;

    // Tracking stats
    size_t max_simul_branches_used = 0;
    size_t total_branches = 0;
    size_t iterations = 0;

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

    // Perform the extension. The while conditions are heuristics to avoid searching
    // the graph too much 
    while(!queue.empty() && iterations++ < MAX_ITERATIONS && queue.size() < MAX_SIMULTANEOUS_BRANCHES && total_branches < MAX_TOTAL_BRANCHES)
    {
        if(queue.size() > max_simul_branches_used)
            max_simul_branches_used = queue.size();

        BuilderExtensionNode curr = queue.front();
        queue.pop();

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
            
            // Check if this sequence is present in the FM-index of the target
            // If so, it is the join point of the de Bruijn graph and we extend no further.
            size_t targetCount = BWTAlgorithms::countSequenceOccurrences(newStr, m_parameters.baseIndex);

            if(targetCount >= MIN_TARGET_COUNT)
            {
                if(curr.direction == ED_SENSE)
                    sense_join_vector.push_back(pVertex);
                else
                    antisense_join_vector.push_back(pVertex);
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
    
    delete pGraph;
    return HBRC_OK;
}
