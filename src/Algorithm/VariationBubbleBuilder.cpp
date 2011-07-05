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

//
//
//
VariationBubbleBuilder::VariationBubbleBuilder() : m_kmerThreshold(1)
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

// The source string is the string the bubble starts from
void VariationBubbleBuilder::setSourceString(const std::string& str)
{
    // Create a new vertex for the source sequence
    // As we are creating a de Bruijn graph, we use the sequence
    // of the vertex as its ID
    Vertex* pVertex = new(m_pGraph->getVertexAllocator()) Vertex(str, str);
    pVertex->setColor(SOURCE_COLOR);
    m_pGraph->addVertex(pVertex);

    // Add the vertex to the extension queue
    m_queue.push(BubbleExtensionNode(pVertex, ED_SENSE));
    m_queue.push(BubbleExtensionNode(pVertex, ED_ANTISENSE));
}

// The source index is the index that the contains the source string
void VariationBubbleBuilder::setSourceIndex(const BWT* pBWT, const BWT* pRBWT)
{
    m_pSourceBWT = pBWT;
    m_pSourceRevBWT = pRBWT;
}

// The target index is the index that we try to build the bubble onto
void VariationBubbleBuilder::setTargetIndex(const BWT* pBWT, const BWT* pRBWT)
{
    m_pTargetBWT = pBWT;
    m_pTargetRevBWT = pRBWT;
}

// Run the bubble construction process
BubbleResult VariationBubbleBuilder::run()
{
    BubbleResult result;
    result.returnCode = BRC_UNKNOWN;

    // Build the source half of the bubble
    result.returnCode = buildSourceBubble();
    if(result.returnCode != BRC_OK)
        return result;
    
    // Build the target half of the bubble
    result.returnCode = buildTargetBubble();
    if(result.returnCode != BRC_OK)
        return result;

    parseBubble(result);
    return result;
}

// Build the portion of the graph from the source vertex
// until it meets the target graph. Returns true
// if a path to the target graph was found.
BubbleResultCode VariationBubbleBuilder::buildSourceBubble()
{
    assert(!m_queue.empty());
    while(!m_queue.empty())
    {
        BubbleExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        std::string extensions = BWTAlgorithms::calculateDeBruijnExtensions(vertStr, m_pSourceBWT, m_pSourceRevBWT, curr.direction, m_kmerThreshold);

        if(extensions.size() > 1)
            return BRC_SOURCE_BRANCH;

        for(size_t i = 0; i < extensions.size(); ++i)
        {
            std::string newStr = makeDeBruijnVertex(vertStr, extensions[i], curr.direction);
            
            // Create the new vertex and edge in the graph

            // If this vertex already exists, the graph must contain a loop
            if(m_pGraph->getVertex(newStr) != NULL)
                return BRC_SOURCE_BRANCH;

            Vertex* pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
            pVertex->setColor(SOURCE_COLOR);
            m_pGraph->addVertex(pVertex);
            addDeBruijnEdges(curr.pVertex, pVertex, curr.direction);
            
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
            }
            else
            {
                // Add the vertex to the extension queue
                m_queue.push(BubbleExtensionNode(pVertex, curr.direction));
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
    m_queue.push(BubbleExtensionNode(m_antisenseJoins.front(), ED_SENSE));

    while(!m_queue.empty())
    {
        BubbleExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        std::string extensions = BWTAlgorithms::calculateDeBruijnExtensions(vertStr, m_pTargetBWT, m_pTargetRevBWT, curr.direction, m_kmerThreshold);

        if(extensions.size() > 1)
            return BRC_TARGET_BRANCH;

        for(size_t i = 0; i < extensions.size(); ++i)
        {
            std::string newStr = makeDeBruijnVertex(vertStr, extensions[i], curr.direction);
            Vertex* pVertex = m_pGraph->getVertex(newStr);
            bool joinFound = false;
            if(pVertex == NULL)
            {
                // Not a join vertex, create a new vertex and add it to the graph and queue
                
                // If this vertex already exists, the graph must contain a loop
                if(m_pGraph->getVertex(newStr) != NULL)
                    return BRC_TARGET_BRANCH;

                pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
                pVertex->setColor(TARGET_COLOR);
                m_pGraph->addVertex(pVertex);
                // Add the vertex to the extension queue
                m_queue.push(BubbleExtensionNode(pVertex, curr.direction));
            }
            else
            {
                if(pVertex->getColor() != JOIN_COLOR)
                {
                    // Vertex exists but it is not the join vertex
                    // This means a simple loop has been found in the target
                    return BRC_TARGET_BRANCH;
                }

                joinFound = true;
            }
            
            // Create the new edge in the graph        
            addDeBruijnEdges(curr.pVertex, pVertex, curr.direction);
            if(joinFound)
                return BRC_OK;
        }
    }

    // no path found
    return BRC_TARGET_BROKEN;
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
                                       10000000, // max distance to search
                                       10000000, // max nodes to search
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
    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        std::string walkStr = outWalks[i].getString(SGWT_START_TO_END);
        bool isTarget = classifyWalk(outWalks[i]);
        if(isTarget)
            targetStrings.push_back(walkStr);
        else
            sourceStrings.push_back(walkStr);
    }
    
    if(targetStrings.size() == 1 && sourceStrings.size() == 1)
    {   
        result.returnCode = BRC_OK;
        result.targetString = targetStrings.front();
        result.sourceString = sourceStrings.front();
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
bool VariationBubbleBuilder::classifyWalk(const SGWalk& walk) const
{
    GraphColor branchCol = GC_WHITE;

    size_t numVertices = walk.getNumVertices();
    if(numVertices <= 2)
    {
        std::cerr << "VariationBubbleBuilder error: degenerate bubble found\n";
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < numVertices; ++i)
    {
        const Vertex* pVertex = walk.getVertex(i);
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
        {
            branchCol = vertCol;
        }
    }

    return branchCol == TARGET_COLOR;
}

// Add a new edge to the graph denoting the relationship between pX and pY.
// Assumes pX and pY are already present in the m_pGraph
void VariationBubbleBuilder::addDeBruijnEdges(const Vertex* pX, const Vertex* pY, EdgeDir direction)
{
    assert(pX->getSeq().length() == pY->getSeq().length());
    
    // overlap length for a de bruijn edge
    size_t p = pX->getSeq().length() - 1;

    // Construct an overlap object for this relationship
    Overlap o;
    o.id[0] = pX->getID();
    o.id[1] = pY->getID();

    o.match.isReverse = false;
    o.match.numDiff = 0;

    if(direction == ED_SENSE)
    {
        // pX -> pY
        o.match.coord[0].interval.start = 1;
        o.match.coord[1].interval.start = 0;
    }
    else
    {
        // pY -> pX
        o.match.coord[0].interval.start = 0;
        o.match.coord[1].interval.start = 1;
    }

    o.match.coord[0].interval.end = o.match.coord[0].interval.start + p - 1; // inclusive coordinate
    o.match.coord[1].interval.end = o.match.coord[1].interval.start + p - 1;
    o.match.coord[0].seqlen = p + 1;
    o.match.coord[1].seqlen = p + 1;
    Edge* e = SGAlgorithms::createEdgesFromOverlap(m_pGraph, o, false);
    assert(e != NULL);
}


// Make the sequence of a new deBruijn vertex using the edge details
std::string VariationBubbleBuilder::makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction)
{
    std::string w;
    size_t p = v.size() - 1;
    if(direction == ED_SENSE)
    {
        w = v.substr(1, p);
        w.append(1, edgeBase);
    }
    else
    {
        w.append(1, edgeBase);
        w.append(v.substr(0, p));
    }
    return w;
}

