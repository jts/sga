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

    // We search until we find the first common vertex in each direction
    size_t MIN_TARGET_COUNT = m_parameters.bReferenceMode ? 1 : 2;
    size_t MAX_ITERATIONS = 2000;
    size_t MAX_SIMULTANEOUS_BRANCHES = 40;
    size_t MAX_TOTAL_BRANCHES = 50;
    int MAX_DISTANCE = 100;

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
        AlphaCount64 extensionCounts = 
            BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_parameters.variantIndex.pBWT, curr.direction);

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

            // Add the vertex to the extension queue
            BuilderExtensionNode nextNode(pVertex, curr.direction, curr.distance + 1);
            if(nextNode.distance < MAX_DISTANCE)
                queue.push(nextNode);
        }
        
        // Update the total number of times we branches the search
        if(!extensionsUsed.empty())
            total_branches += extensionsUsed.size() - 1;
    }

    printf("AS: %zu S: %zu\n", antisense_join_vector.size(), sense_join_vector.size());

    removeNonJoinTips(pGraph, antisense_join_vector, sense_join_vector);
 
    cullRedundantJoinVertices(antisense_join_vector, pGraph, !ED_ANTISENSE);
    cullRedundantJoinVertices(sense_join_vector, pGraph, !ED_SENSE);

    printf("AS: %zu S: %zu\n", antisense_join_vector.size(), sense_join_vector.size());

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

// This algorithm removes tip branches that consist of non-join (tumor-only) vertices
void DeBruijnHaplotypeBuilder::removeNonJoinTips(StringGraph* pGraph, std::vector<Vertex*>& antisense_joins, std::vector<Vertex*>& sense_joins)
{
    pGraph->setColors(GC_WHITE);

    // Mark all join vertices as RED. These are ignored in the subsequent algorithm
    for(size_t i = 0; i < antisense_joins.size(); ++i)
        antisense_joins[i]->setColor(GC_RED);
    for(size_t i = 0; i < sense_joins.size(); ++i)
        sense_joins[i]->setColor(GC_RED);

    VertexPtrVec vertices = pGraph->getAllVertices();
    for(size_t i = 0; i < vertices.size(); ++i)
    {
        Vertex* x = vertices[i];
        _checkNonJoinForTrim(x, ED_SENSE);
        _checkNonJoinForTrim(x, ED_ANTISENSE);
    }

    pGraph->sweepVertices(GC_BLACK);
    pGraph->setColors(GC_BLACK);
}

// Recursive call to trim non-join tip vertices from the graph
void DeBruijnHaplotypeBuilder::_checkNonJoinForTrim(Vertex* x, EdgeDir direction)
{
    // join vertices (red) and previous visited (black) are ignored
    if(x->getColor() == GC_RED || x->getColor() == GC_BLACK)
        return;

    // Check whether there are any edges in this direction
    // Vertices that are marked for deletion (BLACK) are ignored
    EdgePtrVec edges = x->getEdges(direction);

    size_t n_edges = edges.size();
    size_t n_trimmed_edges = 0;
    for(size_t i = 0; i < edges.size(); ++i)
        n_trimmed_edges += (edges[i]->getEnd()->getColor() == GC_BLACK);

    if(n_edges == n_trimmed_edges)
    {
        // this vertex is a tip, remove it
        x->setColor(GC_BLACK);
        printf("Marking for deletion\n");
        // Recurse to neighbors
        for(size_t i = 0; i < edges.size(); ++i)
            _checkNonJoinForTrim(edges[i]->getEnd(), direction);
    }
}

void DeBruijnHaplotypeBuilder::cullRedundantJoinVertices(std::vector<Vertex*>& join_vertices, StringGraph* pGraph, EdgeDir direction)
{
    pGraph->checkColors(GC_BLACK);
    // Remove join nodes that lie on the unipath of each join vertex
    for(size_t i = 0; i < join_vertices.size(); ++i)
    {
        std::vector<Vertex*> unipath(1, join_vertices[i]);

        // Skip vertices already marked
        if(unipath.back()->getColor() == GC_RED)
            continue;

        // Extend unipath
        while(1)
        {
            EdgePtrVec edges = unipath.back()->getEdges(direction);
            if(edges.size() != 1)
                break;
            unipath.push_back(edges[0]->getEnd());
        }

        // Mark the non-terminal vertices of the unipath
        for(size_t j = 1; j < unipath.size(); ++j)
            unipath[j]->setColor(GC_RED);
    }

    std::vector<Vertex*> trimmed;
    for(size_t i = 0; i < join_vertices.size(); ++i)
    {
        if(join_vertices[i]->getColor() != GC_RED)
            trimmed.push_back(join_vertices[i]);
        join_vertices[i]->setColor(GC_BLACK);
    }

    pGraph->setColors(GC_BLACK);
    join_vertices.swap(trimmed);
}

Vertex* DeBruijnHaplotypeBuilder::extendUnambiguously(Vertex* pStart, StringGraph* pGraph, 
                                                      EdgeDir direction, size_t max_distance, size_t min_target_count)
{
    size_t distance = 0;
    Vertex* pCurrent = pStart;

    while(distance < max_distance)
    {
        std::string vertStr = pCurrent->getSeq().toString();
        AlphaCount64 extensionCounts = 
            BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_parameters.variantIndex.pBWT, direction);
        
        // Count valid extensions
        char ext_base = '\0';
        size_t ext_count = 0;
        size_t n_ext = 0;
        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            if(count >= m_parameters.minDBGCount)
            {
                ext_base = b;
                ext_count = count;
                n_ext++;
            }
        }
        
        // Stop if no extension or ambiguous
        if(n_ext != 1)
            return pCurrent;

        std::string newStr = VariationBuilderCommon::makeDeBruijnVertex(vertStr, ext_base, direction);
        
        // Stop if the vertex is already in the graph
        if(pGraph->getVertex(newStr) != NULL)
            return pCurrent;

        // Stop if the vertex is not reprsented in the base sequence
        size_t targetCount = BWTAlgorithms::countSequenceOccurrences(newStr, m_parameters.baseIndex);
        if(targetCount < min_target_count)
            return pCurrent;

        // Create the new vertex and edge in the graph
        // Allocate the new vertex and add it to the graph
        Vertex* pNew = new(pGraph->getVertexAllocator()) Vertex(newStr, newStr);
        pNew->setColor(GC_BLACK);
        pGraph->addVertex(pNew);

        // Add edges
        VariationBuilderCommon::addSameStrandDeBruijnEdges(pGraph, pCurrent, pNew, direction);

        pCurrent = pNew;
        distance++;
    }

    return pCurrent;
}
