///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphCompare - Compare two (abstractly represented)
// graphs against each other to find strings
// only present in one.
//
// The graphs are abstractly represented as
// an FM-index.
//
#include "GraphCompare.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"

//
//
//
void GraphCompareStackNode::initialize(char b, const BWTVector& bwts, const BWTVector& rbwts)
{
    assert(bwts.size() == NUM_GRAPHS && rbwts.size() == NUM_GRAPHS);

    // Initialize the intervals for each BWT to the range containing all
    // suffixes ending with b
    for(size_t i = 0; i < NUM_GRAPHS; ++i)
    {
        BWTAlgorithms::initIntervalPair(intervalPairs[i], b, bwts[i], rbwts[i]);
        assert(intervalPairs[i].isValid());

        // Initialize the extension counts
        lowerCounts[i] = bwts[i]->getFullOcc(intervalPairs[i].interval[0].lower - 1);
        upperCounts[i] = bwts[i]->getFullOcc(intervalPairs[i].interval[0].upper);
    }

    str.append(1,b);
    length = 1;
    alphaIndex = 0;
}

// Update the intervals
void GraphCompareStackNode::update(char b, const BWTVector& bwts, const BWTVector& rbwts)
{
    assert(bwts.size() == NUM_GRAPHS && rbwts.size() == NUM_GRAPHS);

    // Update each interval for the extension symbol b
    for(size_t i = 0; i < NUM_GRAPHS; ++i)
    {
        // Update the interval if it is valid
        if(intervalPairs[i].interval[0].isValid())
            BWTAlgorithms::updateBothL(intervalPairs[i], b, bwts[i], lowerCounts[i], upperCounts[i]);
        
        // Update occurrence counts
        if(intervalPairs[i].interval[0].isValid())
        {
            lowerCounts[i] = bwts[i]->getFullOcc(intervalPairs[i].interval[0].lower - 1);
            upperCounts[i] = bwts[i]->getFullOcc(intervalPairs[i].interval[0].upper);    
        }
        else
        {
            lowerCounts[i].clear();
            upperCounts[i].clear();
        }
    }

    str.append(1,b);
    length += 1;
    alphaIndex = 0;
}

//
AlphaCount64 GraphCompareStackNode::getAggregateExtCount() const
{
    AlphaCount64 out;
    for(size_t i = 0; i < NUM_GRAPHS; ++i)
    {   
        if(intervalPairs[i].interval[0].isValid())
            out += upperCounts[i] - lowerCounts[i];
    }
    return out;
}

//
void GraphCompareStackNode::print() const
{
    std::cout << "Stack string: " << str << "\n";
}


//
//
//
GraphCompare::GraphCompare(const BWT* pBaseBWT, 
                           const BWT* pBaseRBWT,
                           const BWT* pVariantBWT, 
                           const BWT* pVariantRBWT, 
                           int kmer) : m_pBaseBWT(pBaseBWT),  
                                       m_pBaseRevBWT(pBaseRBWT),
                                       m_pVariantBWT(pVariantBWT),  
                                       m_pVariantRevBWT(pVariantRBWT),
                                       m_kmer(kmer)
{
    m_pUsedVariantKmers = new BitVector(pVariantBWT->getBWLen());
}

//
GraphCompare::~GraphCompare()
{
    delete m_pUsedVariantKmers;
}

// Run the actual comparison
void GraphCompare::run()
{
    std::cout << "Running graph comparison\n";

    // Make a vector of the forward and reverse BWTS
    BWTVector bwts;
    bwts.push_back(m_pBaseBWT);
    bwts.push_back(m_pVariantBWT);

    BWTVector rbwts;
    rbwts.push_back(m_pBaseRevBWT);
    rbwts.push_back(m_pVariantRevBWT);

    // 
    GraphCompareStack stack;

    // Initialize the stack with a node for each of ACGT
    for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
    {
        GraphCompareStackNode node;
        node.initialize(DNA_ALPHABET::getBase(i), bwts, rbwts);
        stack.push(node);
    }

    // Loop over the stack and find variant kmers
    size_t loops = 0;
    size_t maxStack = 0;
    size_t kmersFound = 0;
    size_t variantKmersFound = 0;
    size_t unexpectedKmersFound = 0;
    size_t numPush = 0;
    size_t numPop = 0;

    while(!stack.empty())
    {

        GraphCompareStackNode* pNode = &stack.top();
        
        // Progress counter bookkeeping
        loops += 1;
        if(stack.size() > maxStack)
            maxStack = stack.size();

        if(loops % 10000000 == 0)
        {
            std::cout << "Loop: " << loops << "\n";
            std::cout << "Kmers: " << kmersFound << "\n";
            //pNode->print();
        }

        bool doPop = true;

        if(pNode->length == (int)m_kmer)
        {
            // do something
            kmersFound += 1;

            if(!pNode->intervalPairs[0].isValid() && pNode->intervalPairs[1].isValid())
            {
                if(processVariantKmer(reverse(pNode->str), bwts, rbwts, 1))
                    variantKmersFound += 1;
            }
            
            /*
            if(pNode->intervalPairs[0].isValid() && !pNode->intervalPairs[1].isValid())
            {
                if(processVariantKmer(reverse(pNode->str), bwts, rbwts))
                    unexpectedKmersFound += 1;
            }
            */
        }
        else
        {
            // extend
            // Calculate the aggregate extension count
            AlphaCount64 ext_count = pNode->getAggregateExtCount(); 
            bool hasBranch = ext_count.getNumNonZero() > 0 && !ext_count.hasUniqueDNAChar();

            while(pNode->alphaIndex < DNA_ALPHABET::size)
            {
                char b = DNA_ALPHABET::getBase(pNode->alphaIndex);
                if(ext_count.get(b) > 0)
                {
                    // Branch the search by creating a new stack entry
                    if(hasBranch)
                    {
                        GraphCompareStackNode branched = *pNode;
                        branched.update(b, bwts, rbwts);
                        pNode->alphaIndex += 1;
                        stack.push(branched);
                        numPush += 1;
                        doPop = false;
                    }
                    else
                    {
                        // Update the interval for the current stack entry
                        pNode->update(b, bwts, rbwts);
                        doPop = false;
                    }
                    break;
                }
                else
                {
                    pNode->alphaIndex += 1;
                }
            }
        }

        if(doPop)
        {
            pNode = NULL;
            stack.pop();
            numPop += 1;
        }
    }

    printf("Done traversal: %zu\n", loops);
    printf("Max stack: %zu\n", maxStack);
    printf("Push: %zu\n", numPush);
    printf("Pop: %zu\n", numPop);
    printf("Total kmers: %zu\n", kmersFound);
    printf("Variant kmers: %zu\n", variantKmersFound);
    printf("Unexpected kmers: %zu\n", unexpectedKmersFound);

}   

//
bool GraphCompare::processVariantKmer(const std::string& str, const BWTVector& bwts, const BWTVector& rbwts, int varIndex)
{
    assert(varIndex == 0 || varIndex == 1);
    // Check the count of the kmer in both intervals
    size_t c1 = BWTAlgorithms::countSequenceOccurrences(str, bwts[0]);
    size_t c2 = BWTAlgorithms::countSequenceOccurrences(str, bwts[1]);
    
    // If there is reverse-complement coverage of this kmer, do not use
    // it as a bubble seed
    if(c1 > 0 && c2 > 0)
        return false;

    // Check if this k-mer has been marked in the bitvector
    bool seen = isKmerMarked(str);
    if(seen)
        return false;

    BubbleBuilder builder;
    builder.setSourceIndex(bwts[varIndex], rbwts[varIndex]);
    builder.setTargetIndex(bwts[1 - varIndex], rbwts[1 - varIndex]);
    builder.setSourceString(str);

    //
    BubbleResult result = builder.run();
    if(result.success)
    {
        assert(!result.targetString.empty());
        assert(!result.sourceString.empty());
        StdAlnTools::globalAlignment(result.targetString, result.sourceString, true);
        markVariantSequenceKmers(result.sourceString);
    }
    
    return result.success;
}

// Update the bit vector with the kmers that were assembled into str
void GraphCompare::markVariantSequenceKmers(const std::string& str)
{
    assert(str.size() >= m_kmer);
    size_t n = str.size() - m_kmer + 1;

    for(size_t i = 0; i < n; ++i)
    {
        std::string kseq = str.substr(i, m_kmer);
        BWTInterval interval = BWTAlgorithms::findInterval(m_pVariantBWT, kseq);
        if(interval.isValid())
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
                m_pUsedVariantKmers->updateCAS(j, false, true);
        }

        // Mark the reverse complement k-mers too
        std::string rc_kseq = reverseComplement(kseq);
        interval = BWTAlgorithms::findInterval(m_pVariantBWT, rc_kseq);
        if(interval.isValid())
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
                m_pUsedVariantKmers->updateCAS(j, false, true);
        }
    }
}

//
bool GraphCompare::isKmerMarked(const std::string& str) const
{
    assert(str.size() == m_kmer);
    BWTInterval interval = BWTAlgorithms::findInterval(m_pVariantBWT, str);
    if(interval.isValid())
    {
        if(m_pUsedVariantKmers->test(interval.lower))
            return true;
    }
    
    // Mark the reverse complement k-mers too
    std::string rc_str = reverseComplement(str);
    interval = BWTAlgorithms::findInterval(m_pVariantBWT, rc_str);
    if(interval.isValid())
    {
        if(m_pUsedVariantKmers->test(interval.lower))
            return true;
    }
    return false;
}

//
//
//
BubbleBuilder::BubbleBuilder()
{
    m_pGraph = new StringGraph;
}

//
BubbleBuilder::~BubbleBuilder()
{
    delete m_pGraph;
}

// The source string is the string the bubble starts from
void BubbleBuilder::setSourceString(const std::string& str)
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
void BubbleBuilder::setSourceIndex(const BWT* pBWT, const BWT* pRBWT)
{
    m_pSourceBWT = pBWT;
    m_pSourceRevBWT = pRBWT;
}

// The target index is the index that we try to build the bubble onto
void BubbleBuilder::setTargetIndex(const BWT* pBWT, const BWT* pRBWT)
{
    m_pTargetBWT = pBWT;
    m_pTargetRevBWT = pRBWT;
}

// Run the bubble construction process
BubbleResult BubbleBuilder::run()
{
    std::cout << "Running bubble builder\n";
    BubbleResult result;
    result.success = false;

    // Build the source half of the bubble
    result.success = buildSourceBubble();
    if(!result.success)
        return result;
    
    // Build the target half of the bubble
    result.success = buildTargetBubble();
    if(!result.success)
    {
        std::cout << "Failed to build target\n";
        return result;
    }

    result = parseBubble();
    return result;
}

// Build the portion of the graph from the source vertex
// until it meets the target graph. Returns true
// if a path to the target graph was found.
bool BubbleBuilder::buildSourceBubble()
{
    assert(!m_queue.empty());
    while(!m_queue.empty())
    {
        BubbleExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        std::string extensions = BWTAlgorithms::calculateDeBruijnExtensions(vertStr, m_pSourceBWT, m_pSourceRevBWT, curr.direction);

        if(extensions.size() > 1)
            std::cout << "Branching\n";

        for(size_t i = 0; i < extensions.size(); ++i)
        {
            std::string newStr = makeDeBruijnVertex(vertStr, extensions[i], curr.direction);
            
            // Create the new vertex and edge in the graph
            Vertex* pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
            pVertex->setColor(SOURCE_COLOR);
            m_pGraph->addVertex(pVertex);
            addDeBruijnEdges(curr.pVertex, pVertex, curr.direction);
            
            // Check if this sequence is present in the FM-index of the target
            // If so, it is the join point of the de Bruijn graph and we extend no further.
            size_t targetCount = BWTAlgorithms::countSequenceOccurrences(newStr, m_pTargetBWT);
            if(targetCount > 0)
            {
                std::cout << "Join point found for direction: " << curr.direction << " " << newStr << "\n";
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
        return true;
    else
        return false;
}

// Build the portion of the graph between the found target
// join vertices.
bool BubbleBuilder::buildTargetBubble()
{
    assert(m_queue.empty());
    assert(m_antisenseJoins.size() == 1);
    assert(m_senseJoins.size() == 1);
    std::cout << "TARGET BUBBLE\n";
    // Add the antisense join vertex to the extension queue
    m_queue.push(BubbleExtensionNode(m_antisenseJoins.front(), ED_SENSE));

    while(!m_queue.empty())
    {
        BubbleExtensionNode curr = m_queue.front();
        m_queue.pop();

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        std::string extensions = BWTAlgorithms::calculateDeBruijnExtensions(vertStr, m_pTargetBWT, m_pTargetRevBWT, curr.direction);

        if(extensions.size() > 1)
            return false;

        for(size_t i = 0; i < extensions.size(); ++i)
        {
            std::string newStr = makeDeBruijnVertex(vertStr, extensions[i], curr.direction);
            Vertex* pVertex = m_pGraph->getVertex(newStr);
            bool joinFound = false;
            if(pVertex == NULL)
            {
                // Not a join vertex, create a new vertex and add it to the graph and queue
                pVertex = new(m_pGraph->getVertexAllocator()) Vertex(newStr, newStr);
                pVertex->setColor(TARGET_COLOR);
                m_pGraph->addVertex(pVertex);
                // Add the vertex to the extension queue
                m_queue.push(BubbleExtensionNode(pVertex, curr.direction));
            }
            else
            {
                assert(pVertex->getColor() == JOIN_COLOR);
                joinFound = true;
            }
            
            // Create the new edge in the graph        
            addDeBruijnEdges(curr.pVertex, pVertex, curr.direction);
            if(joinFound)
                return true;
        }
    }

    // no path found
    return false;
}

// After the bubble has been built into the graph, this function
// finds and compares the two sequences
BubbleResult BubbleBuilder::parseBubble()
{
    BubbleResult result;
    result.success = false;
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
        return result;

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
        result.success = true;
        result.targetString = targetStrings.front();
        result.sourceString = sourceStrings.front();
    }
    else
    {
        result.success = false;
    }
    return result;
}

// Returns true if the walk is the part of the target sequence
bool BubbleBuilder::classifyWalk(const SGWalk& walk) const
{
    GraphColor branchCol = GC_WHITE;

    size_t numVertices = walk.getNumVertices();
    if(numVertices <= 2)
    {
        std::cerr << "BubbleBuilder error: degenerate bubble found\n";
        exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < numVertices; ++i)
    {
        const Vertex* pVertex = walk.getVertex(i);
        GraphColor vertCol = pVertex->getColor();
        if(vertCol == JOIN_COLOR && i != 0 && i != numVertices - 1)
        {
            std::cerr << "BubbleBuilder error: interior join vertex found\n";
            exit(EXIT_FAILURE);
        }
        
        if((vertCol == TARGET_COLOR || vertCol == SOURCE_COLOR) && branchCol != GC_WHITE && branchCol != vertCol)
        {
            std::cerr << "BubbleBuilder error: unexpected mixed-color branch\n";
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
void BubbleBuilder::addDeBruijnEdges(const Vertex* pX, const Vertex* pY, EdgeDir direction)
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
std::string BubbleBuilder::makeDeBruijnVertex(const std::string& v, char edgeBase, EdgeDir direction)
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

