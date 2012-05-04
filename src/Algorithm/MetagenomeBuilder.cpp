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

//
MetagenomeBuilder::MetagenomeBuilder()
{
    m_frequencyFilter = 0.5;
    m_hardMinCoverage = 3;
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

    while(!m_queue.empty())
    {
        numIters += 1;
        BuilderExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string strX = curr.pVertex->getSeq().toString();

        // Count the number of branches from this sequence
        std::pair<std::string, int> nodeY = getBestEdgeNode(strX, m_vertexCoverageMap[strX], curr.direction);
        std::string& strY = nodeY.first;
        int coverageY = nodeY.second;

        if(strY != "")
        {
            // Get the best node in the opposite direction for 
            std::pair<std::string, int> nodeZ = getBestEdgeNode(strY, coverageY, !curr.direction);
            std::string& strZ = nodeZ.first;
            if(strZ == "" || strZ != strX)
            {
                // Either Y does not have an unambiguous best connection or it is
                // not X. Either way, we stop the extension.
                continue;
            }
        }
        else
        {
            continue;
        }

        // Create the new vertex and edge in the graph
        // If this vertex already exists, the graph must contain a loop so we stop
        if(m_pGraph->getVertex(strY) != NULL)
            break;

        Vertex* pNewVertex = new(m_pGraph->getVertexAllocator()) Vertex(strY, strY);
        addVertex(pNewVertex, coverageY);
        VariationBuilderCommon::addSameStrandDeBruijnEdges(m_pGraph, curr.pVertex, pNewVertex, curr.direction);
            
        // Add the vertex to the extension queue
        m_queue.push(BuilderExtensionNode(pNewVertex, curr.direction));
    }
    // Done extension
}

// Get the best de Bruijn graph node connected to nodeX
// If ambiguous, the empty string is returned
std::pair<std::string, int> MetagenomeBuilder::getBestEdgeNode(const std::string& nodeX, size_t nodeCoverage, EdgeDir direction)
{
    AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensions(nodeX, 
                                                                              m_pBWT, 
                                                                              m_pRevBWT, 
                                                                              direction,
                                                                              m_pBWTCache, 
                                                                              m_pRevBWTCache);

    size_t cov_threshold = static_cast<size_t>(std::max(m_frequencyFilter * nodeCoverage, (double)m_hardMinCoverage));
    bool uniqueExtension = extensionCounts.hasUniqueDNAChar() || 
                                VariationBuilderCommon::countValidExtensions(extensionCounts, cov_threshold) == 1;

    std::pair<std::string, int> ret;
    // Fail due to ambiguity
    if(!uniqueExtension)
        return ret;

    // Output the sequence of the vertex linked
    for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
    {
        char b = DNA_ALPHABET::getBase(i);
        size_t count = extensionCounts.get(b);
        if(!uniqueExtension || count < cov_threshold)
            continue; 
        ret.first = VariationBuilderCommon::makeDeBruijnVertex(nodeX, b, direction);
        ret.second = count;
        return ret;
    }
    return ret;
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
