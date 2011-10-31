///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// VariationBubbleBuilder - Construct a variation
// bubble from an initial seed k-mer which only
// appears in one out of a pair of abstract
// deBruijn graphs.
//
#include "VariationBubbleBuilder.h"
#include "BWTAlgorithms.h"
#include "SGSearch.h"
#include "SGAlgorithms.h"
#include "BuilderCommon.h"
#include "HaplotypeBuilder.h"

//
//
//
VariationBubbleBuilder::VariationBubbleBuilder() : m_kmerThreshold(1), m_allowedTargetBranches(0)
{
    m_pGraph = new StringGraph;
}

//
VariationBubbleBuilder::~VariationBubbleBuilder()
{
    delete m_pGraph;
}

//
void VariationBubbleBuilder::setKmerThreshold(size_t t)
{
    m_kmerThreshold = t;
}

//
void VariationBubbleBuilder::setAllowedBranches(size_t b)
{
    m_allowedTargetBranches = b;
}

// The source string is the string the bubble starts from
void VariationBubbleBuilder::setSourceString(const std::string& str, int coverage)
{
    // Create a new vertex for the source sequence
    // As we are creating a de Bruijn graph, we use the sequence
    // of the vertex as its ID
    Vertex* pVertex = new(m_pGraph->getVertexAllocator()) Vertex(str, str);
    pVertex->setColor(SOURCE_COLOR);
    addVertex(pVertex, coverage);

    // Add the vertex to the extension queue
    m_queue.push(BuilderExtensionNode(pVertex, ED_SENSE));
    m_queue.push(BuilderExtensionNode(pVertex, ED_ANTISENSE));
}

// The source index is the index that the contains the source string
void VariationBubbleBuilder::setSourceIndex(const BWT* pBWT)
{
    m_pSourceBWT = pBWT;
}

// The target index is the index that we try to build the bubble onto
void VariationBubbleBuilder::setTargetIndex(const BWT* pBWT)
{
    m_pTargetBWT = pBWT;
}

// Run the bubble construction process
BubbleResult VariationBubbleBuilder::run()
{


    BubbleResult result;
    result.returnCode = BRC_UNKNOWN;

    bool useHaplotypeBuilder = false;


    if(!useHaplotypeBuilder)
    {
        // Build the source half of the bubble
        result.returnCode = buildSourceBubble();
        if(result.returnCode != BRC_OK)
            return result;
        
        // Build the target half of the bubble
        result.returnCode = buildTargetBubble();
        if(result.returnCode != BRC_OK)
            return result;

        parseBubble(result);
    }
    else
    {
        // Build the source half of the bubble
        result.returnCode = buildSourceBubble();
        if(result.returnCode != BRC_OK)
            return result;
        
        // Build the target half of the bubble
        result.returnCode = buildTargetBubbleHB();
        if(result.returnCode != BRC_OK)
            return result;

        parseStringsHB(result);
    }
    return result;
}

// Perform a breadth-first search from the variant kmer
// until it meets the common portion of the graph. This 
// function will generate de bruijn graph to the  closest
// common vertices.
BubbleResultCode VariationBubbleBuilder::buildSourceBubble()
{
    assert(!m_queue.empty());

    // We search until we find the first common vertex in each direction
    bool joinFound[ED_COUNT];
    joinFound[ED_SENSE] = false;
    joinFound[ED_ANTISENSE] = false;

    while(!m_queue.empty())
    {
        BuilderExtensionNode curr = m_queue.front();
        m_queue.pop();

        // We have found a join in this direction, stop the search
        if(joinFound[curr.direction])
            continue;

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_pSourceBWT, curr.direction);
        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            bool acceptExt = count >= m_kmerThreshold || (count > 0 && extensionCounts.hasUniqueDNAChar());
            if(!acceptExt)
                continue;

            std::string newStr = BuilderCommon::makeDeBruijnVertex(vertStr, b, curr.direction);
            
            // Create the new vertex and edge in the graph
            // Skip if the vertex already exists
            if(m_pGraph->getVertex(newStr) != NULL)
                continue;

            Vertex* pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
            pVertex->setColor(SOURCE_COLOR);
            addVertex(pVertex, count);
            BuilderCommon::addSameStrandDeBruijnEdges(m_pGraph, curr.pVertex, pVertex, curr.direction);
            
            // Check if this sequence is present in the FM-index of the target
            // If so, it is the join point of the de Bruijn graph and we extend no further.
            size_t targetCount = BWTAlgorithms::countSequenceOccurrences(newStr, m_pTargetBWT);
            if(targetCount > 0)
            {
                pVertex->setColor(JOIN_COLOR);
                if(curr.direction == ED_SENSE)
                    m_senseJoins.push_back(pVertex);
                else
                    m_antisenseJoins.push_back(pVertex);

                joinFound[curr.direction] = true;
            }
            else
            {
                // Add the vertex to the extension queue
                m_queue.push(BuilderExtensionNode(pVertex, curr.direction));
            }
        }
    }

    // Check if a unique join path was found
    if(m_antisenseJoins.size() == 1 && m_senseJoins.size() == 1)
        return BRC_OK;
    else
        return BRC_SOURCE_BROKEN;
}

// Build the portion of the graph between the found target
// join vertices.
BubbleResultCode VariationBubbleBuilder::buildTargetBubble()
{
    assert(m_queue.empty());
    assert(m_antisenseJoins.size() == 1);
    assert(m_senseJoins.size() == 1);

    // Add the antisense join vertex to the extension queue
    m_queue.push(BuilderExtensionNode(m_antisenseJoins.front(), ED_SENSE));

    size_t MAX_SIMULTANEOUS_BRANCHES = 20;

    while(!m_queue.empty())
    {
        if(m_queue.size() > MAX_SIMULTANEOUS_BRANCHES)
            return BRC_TARGET_BRANCH;

        BuilderExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_pTargetBWT, curr.direction);
        
        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            bool acceptExt = count >= m_kmerThreshold || (count > 0 && extensionCounts.hasUniqueDNAChar());
            if(!acceptExt)
                continue;

            std::string newStr = BuilderCommon::makeDeBruijnVertex(vertStr, b, curr.direction);
            Vertex* pVertex = m_pGraph->getVertex(newStr);
            bool joinFound = false;
            if(pVertex == NULL)
            {
                pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
                pVertex->setColor(TARGET_COLOR);
                addVertex(pVertex, count);

                // Add the vertex to the extension queue
                m_queue.push(BuilderExtensionNode(pVertex, curr.direction));
            }
            else
            {
                if(pVertex->getColor() == JOIN_COLOR)
                    joinFound = true;
            }
            
            // Create the new edge in the graph        
            BuilderCommon::addSameStrandDeBruijnEdges(m_pGraph, curr.pVertex, pVertex, curr.direction);

            // If we've found the join vertex, we have completed the target half of the bubble
            if(joinFound)
                return BRC_OK;
        }
    }

    // no path found
    return BRC_TARGET_BROKEN;
}

// Build the target half of the bubble using the haplotype builder
BubbleResultCode VariationBubbleBuilder::buildTargetBubbleHB()
{
    assert(m_queue.empty());
    assert(m_antisenseJoins.size() == 1);
    assert(m_senseJoins.size() == 1);

    // Set the start/end points of the haplotype builder
    std::string startSource = m_antisenseJoins.front()->getSeq().toString();
    std::string endSource = m_senseJoins.front()->getSeq().toString();

    size_t TARGET_KMER = 31;
    assert(startSource.size() >= TARGET_KMER);
    std::string startAnchorSeq = startSource.substr(0, TARGET_KMER);
    std::string endAnchorSeq = endSource.substr(endSource.size() - TARGET_KMER);

    AnchorSequence startAnchor;
    startAnchor.sequence = startAnchorSeq;
    startAnchor.count = 0;
    startAnchor.position = 0;

    AnchorSequence endAnchor;
    endAnchor.sequence = endAnchorSeq;
    endAnchor.count = 0;
    endAnchor.position = 0;

    HaplotypeBuilder builder;
    builder.setTerminals(startAnchor, endAnchor);
    builder.setIndex(m_pTargetBWT, NULL);
    builder.setKmerParameters(TARGET_KMER, m_kmerThreshold);

    // Run the builder
    HaplotypeBuilderReturnCode code = builder.run();
    HaplotypeBuilderResult result;

    // The search was successful, build strings from the walks
    BubbleResultCode bubbleCode = BRC_HB_FAILED;
    if(code == HBRC_OK)
    {
        code = builder.parseWalks(result);
        if(code == HBRC_OK)
        {
            m_targetStrings = result.haplotypes;
            bubbleCode = BRC_OK;
        }
    }

    return bubbleCode;
}

// After the bubble has been built into the graph, this function
// finds and compares the two sequences
void VariationBubbleBuilder::parseBubble(BubbleResult& result)
{
    // Parse walks from the graph that go through the bubbles
    SGWalkVector outWalks;
    bool success = SGSearch::findWalks(m_antisenseJoins.front(),
                                       m_senseJoins.front(),
                                       ED_SENSE,
                                       100000, // max distance to search
                                       10000, // max nodes to search
                                       true, // exhaustive search
                                       outWalks);
    if(!success)
    {
        result.returnCode = BRC_WALK_FAILED;
        return;
    }

    // Convert the walks into strings
    StringVector sourceStrings;
    StringVector targetStrings;
    DoubleVector sourceCoverages;
    DoubleVector targetCoverages;

    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        std::string walkStr = outWalks[i].getString(SGWT_START_TO_END);
        int walkCoverage = 0;
        bool isTarget = classifyWalk(outWalks[i], walkCoverage);
        if(isTarget)
        {
            targetStrings.push_back(walkStr);
            targetCoverages.push_back((double)walkCoverage / outWalks[i].getNumVertices());
        }
        else
        {
            sourceStrings.push_back(walkStr);
            sourceCoverages.push_back((double)walkCoverage / outWalks[i].getNumVertices());
        }
    }
    
    if(targetStrings.size() == 1 && sourceStrings.size() == 1)
    {   
        result.returnCode = BRC_OK;
        result.targetString = targetStrings.front();
        result.sourceString = sourceStrings.front();
        result.targetCoverage = targetCoverages.front();
        result.sourceCoverage = sourceCoverages.front();
    }
    else
    {
        result.returnCode = BRC_NO_SOLUTION;
    }
    return;
}

// Parse strings from the graph and the haplotype builder process
void VariationBubbleBuilder::parseStringsHB(BubbleResult& result)
{
    // Parse the source strings directly out of the graph
    SGWalkVector outWalks;
    bool success = SGSearch::findWalks(m_antisenseJoins.front(),
                                       m_senseJoins.front(),
                                       ED_SENSE,
                                       100000, // max distance to search
                                       10000, // max nodes to search
                                       true, // exhaustive search
                                       outWalks);
    if(!success)
    {
        result.returnCode = BRC_WALK_FAILED;
        return;
    }

    // Convert the walks into strings
    StringVector sourceStrings;
    DoubleVector sourceCoverages;

    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        std::string walkStr = outWalks[i].getString(SGWT_START_TO_END);
        int walkCoverage = 0;
        classifyWalk(outWalks[i], walkCoverage);
        sourceStrings.push_back(walkStr);
        sourceCoverages.push_back((double)walkCoverage / outWalks[i].getNumVertices());
    }
    
    if(sourceStrings.size() == 1 && m_targetStrings.size() == 1)
    {   
        result.returnCode = BRC_OK;
        result.sourceCoverage = sourceCoverages.front();
        result.sourceString = sourceStrings.front();
        result.targetString = m_targetStrings.front();
        result.targetCoverage = -1;
    }
    else
    {
        result.returnCode = BRC_NO_SOLUTION;
    }
    return;
}

// Return a vector of kmers on the source portion of the graph
StringVector VariationBubbleBuilder::getSourceKmers() const
{
    StringVector out;
    VertexPtrVec vertexPtrs =  m_pGraph->getAllVertices();
    for(size_t i = 0; i < vertexPtrs.size(); ++i)
    {
        if(vertexPtrs[i]->getColor() == SOURCE_COLOR)
            out.push_back(vertexPtrs[i]->getSeq().toString());
    }
    return out;
}

// Returns true if the walk is the part of the target sequence
// The total coverage of the walk is written to outCoverage
bool VariationBubbleBuilder::classifyWalk(const SGWalk& walk, int& outCoverage) const
{
    GraphColor branchCol = GC_WHITE;

    size_t numVertices = walk.getNumVertices();
    if(numVertices <= 2)
    {
        std::cerr << "VariationBubbleBuilder warning: degenerate bubble found\n";
        return false;
    }

    outCoverage = 0;
    for(size_t i = 0; i < numVertices; ++i)
    {
        const Vertex* pVertex = walk.getVertex(i);

        // Update color state
        GraphColor vertCol = pVertex->getColor();
        if(vertCol == JOIN_COLOR && i != 0 && i != numVertices - 1)
        {
            std::cerr << "VariationBubbleBuilder error: interior join vertex found\n";
            exit(EXIT_FAILURE);
        }
        
        if((vertCol == TARGET_COLOR || vertCol == SOURCE_COLOR) && branchCol != GC_WHITE && branchCol != vertCol)
        {
            std::cerr << "VariationBubbleBuilder error: unexpected mixed-color branch\n";
            std::cerr << "BranchColor: " << Bigraph::getColorString(branchCol) << " vertColor: " << Bigraph::getColorString(vertCol) << "\n";
            exit(EXIT_FAILURE);
        }
        
        if(vertCol == TARGET_COLOR || vertCol == SOURCE_COLOR)
            branchCol = vertCol;

        // Update coverage
        StrIntMap::const_iterator iter = m_vertexCoverageMap.find(pVertex->getSeq().toString());
        assert(iter != m_vertexCoverageMap.end());
        outCoverage += iter->second;
    }
    return branchCol == TARGET_COLOR;
}

// Add a vertex to the graph
void VariationBubbleBuilder::addVertex(Vertex* pVertex, int coverage)
{
    m_pGraph->addVertex(pVertex);
    m_vertexCoverageMap[pVertex->getSeq().toString()] = coverage;
}
