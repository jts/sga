///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// MetagenomeBuilder -- Implementation of the metagenome
// assembly process. This class takes in a kmer and 
// assembles a contig locally around that sequence
//
#include "MetagenomeBuilder.h"
#include "BWTAlgorithms.h"
#include "BuilderCommon.h"

//
MetagenomeBuilder::MetagenomeBuilder()
{
    m_pGraph = new StringGraph;
}

//
MetagenomeBuilder::~MetagenomeBuilder()
{
    delete m_pGraph;
}

void MetagenomeBuilder::setSource(const std::string& seq, int coverage)
{
    Vertex* pVertex = new(m_pGraph->getVertexAllocator()) Vertex(seq, seq);
    addVertex(pVertex, coverage);

    // Add the vertex to the extension queue
    m_queue.push(BuilderExtensionNode(pVertex, ED_SENSE));
    m_queue.push(BuilderExtensionNode(pVertex, ED_ANTISENSE));
}

//
void MetagenomeBuilder::setIndex(const BWT* pBWT, const BWT* pRevBWT, 
                                 const BWTIntervalCache* pFwdCache, const BWTIntervalCache* pRevCache)
{
    m_pBWT = pBWT;
    m_pRevBWT = pRevBWT;
    m_pBWTCache = pFwdCache;
    m_pRevBWTCache = pRevCache;
}

//
void MetagenomeBuilder::setKmerParameters(size_t k, size_t threshold)
{
    m_kmer = k;
    m_kmerThreshold = threshold;
}

//
void MetagenomeBuilder::run()
{
    assert(!m_queue.empty());
    size_t numIters = 0;
    double frequencyFilter = 0.2;

    while(!m_queue.empty())
    {
        numIters += 1;
        BuilderExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensions(vertStr, 
                                                                                  m_pBWT, 
                                                                                  m_pRevBWT, 
                                                                                  curr.direction,
                                                                                  m_pBWTCache, 
                                                                                  m_pRevBWTCache);

        // Count the number of branches from this sequence
        size_t num_branches = BuilderCommon::filterLowFrequency(extensionCounts, frequencyFilter);

        // Fail due to a high-coverage split occuring
        if(num_branches > 1)
            continue;

        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            if(count == 0)
                continue;

            std::string newStr = BuilderCommon::makeDeBruijnVertex(vertStr, b, curr.direction);
         
            // Check if the new sequence to be added into the graph branches in the opposite
            // direction of the assembly. If so, we are entering a repeat and want to stop
            AlphaCount64 extensionCountsIn = BWTAlgorithms::calculateDeBruijnExtensions(vertStr, 
                                                                                        m_pBWT, 
                                                                                        m_pRevBWT, 
                                                                                        !curr.direction,
                                                                                        m_pBWTCache, 
                                                                                        m_pRevBWTCache);

            size_t num_branches_in = BuilderCommon::filterLowFrequency(extensionCountsIn, frequencyFilter);
            if(num_branches_in > 1)
                continue;

            // Create the new vertex and edge in the graph
            // If this vertex already exists, the graph must contain a loop
            if(m_pGraph->getVertex(newStr) != NULL)
                continue; // This vertex exists, a loop has been found. Stop extension

            Vertex* pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
            addVertex(pVertex, count);
            BuilderCommon::addSameStrandDeBruijnEdges(m_pGraph, curr.pVertex, pVertex, curr.direction);
            
            // Add the vertex to the extension queue
            m_queue.push(BuilderExtensionNode(pVertex, curr.direction));
        }
    }
    // Done extension
}

void MetagenomeBuilder::getContigs(StringVector& contigs)
{
    m_pGraph->simplify();
    m_pGraph->getVertexSequences(contigs);
}


// Add a vertex to the graph
void MetagenomeBuilder::addVertex(Vertex* pVertex, int coverage)
{
    m_pGraph->addVertex(pVertex);
    m_vertexCoverageMap[pVertex->getSeq().toString()] = coverage;
}
