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
// KmerRange
//
void KmerRange::initializeBatch(int batchIdx, int totalBatches)
{
    assert(batchIdx < totalBatches);
    int totalKmers = 1 << 2*keyLength;

    int kmersPerBatch = totalKmers / totalBatches;
    start = batchIdx * kmersPerBatch;
    end = start + kmersPerBatch;
}

// Return the kmer cooresponding to the given value
std::string KmerRange::decode(int64_t i) const
{
    assert(i >= start && i < end);
    std::string out(keyLength, 'A');
    for(int64_t k = 0; k < keyLength; ++k)
    {
        // Get the character as position k (0 = left-most position)
        size_t code = (i >> 2*(keyLength - k - 1)) & 3;
        char b = DNA_ALPHABET::getBase(code);
        out[k] = b;
    }
    return out;
}

//
// GraphCompareStats
//
void GraphCompareStats::clear()
{
    numBubbles = 0;
    numAttempted = 0;
    numTargetBranched = 0;
    numSourceBranched = 0;
    numTargetBroken = 0;
    numSourceBroken = 0;
    numWalkFailed = 0;
    numNoSolution = 0;

    numInsertions = 0;
    numDeletions = 0;
    numSubs = 0;
}

//
void GraphCompareStats::add(const GraphCompareStats& other)
{
    numBubbles += other.numBubbles;
    numAttempted += other.numAttempted;
    numTargetBranched += other.numTargetBranched;
    numSourceBranched += other.numSourceBranched;
    numTargetBroken += other.numTargetBroken;
    numSourceBroken += other.numSourceBroken;
    numWalkFailed += other.numWalkFailed;
    numNoSolution += other.numNoSolution;

    numInsertions += other.numInsertions;
    numDeletions += other.numDeletions;
    numSubs += other.numSubs;
}

//
void GraphCompareStats::print() const
{
    printf("Total attempts: %d\n", numAttempted);
    printf("Total bubbles: %d\n", numBubbles);
    printf("Failed - target branched: %d\n", numTargetBranched);
    printf("Failed - source branched: %d\n", numSourceBranched);
    printf("Failed - target broken: %d\n", numTargetBroken);
    printf("Failed - source broken: %d\n", numSourceBroken);
    printf("Failed - no walk: %d\n", numWalkFailed);
    printf("Failed - no solution: %d\n", numNoSolution);
    
    printf("Num subs found: %d\n", numSubs);
    printf("Num insertions found: %d\n", numInsertions);
    printf("Num deletions found: %d\n", numDeletions);
}

//
// GraphCompareStackNode
//
void GraphCompareStackNode::initialize(const std::string& in, const BWTVector& bwts, const BWTVector& rbwts)
{
    assert(bwts.size() == NUM_GRAPHS && rbwts.size() == NUM_GRAPHS);

    // Initialize the intervals for each BWT to the range containing all
    // strings starting with in
    for(size_t i = 0; i < NUM_GRAPHS; ++i)
    {
        intervalPairs[i] = BWTAlgorithms::findIntervalPair(bwts[i], rbwts[i], in);
        assert(intervalPairs[i].isValid());

        // Initialize the extension counts
        lowerCounts[i] = bwts[i]->getFullOcc(intervalPairs[i].interval[0].lower - 1);
        upperCounts[i] = bwts[i]->getFullOcc(intervalPairs[i].interval[0].upper);
    }

    reverseStr = reverse(in);
    length = in.size();
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

    reverseStr.append(1,b);
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
    std::cout << "Stack string: " << reverse(reverseStr) << "\n";
}

//
//
//
GraphCompare::GraphCompare(const GraphCompareParameters& params) : m_parameters(params)
{
    //m_pWriter = createWriter("variants.fa");
    m_stats.clear();
}

//
GraphCompare::~GraphCompare()
{
    //delete m_pWriter;
}

// Run the actual comparison
void GraphCompare::run()
{
    printf("Running graph comparison for range [%d %d]\n", (int)m_parameters.range.start, (int)m_parameters.range.end);

    // Make a vector of the forward and reverse BWTS
    BWTVector bwts;
    bwts.push_back(m_parameters.pBaseBWT);
    bwts.push_back(m_parameters.pVariantBWT);

    BWTVector rbwts;
    rbwts.push_back(m_parameters.pBaseRevBWT);
    rbwts.push_back(m_parameters.pVariantRevBWT);

    // 
    GraphCompareStack stack;
    initializeStack(m_parameters.range, stack, bwts, rbwts);

    // Loop over the stack and find variant kmers
    size_t loops = 0;
    size_t maxStack = 0;
    size_t kmersFound = 0;
    size_t variantKmersFound = 0;
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
            //std::cout << "Loop: " << loops << "\n";
            //std::cout << "Kmers: " << kmersFound << "\n";
            //pNode->print();
        }


        if(pNode->length == m_parameters.kmer)
        {
            // do something
            kmersFound += 1;
            
            if(isVariantKmer(pNode))
            {
                if(processVariantKmer(reverse(pNode->reverseStr), bwts, rbwts, 1))
                    variantKmersFound += 1;
            }
        }

        bool doPop = updateNodeAndStack(pNode, stack, bwts, rbwts);
        if(doPop)
        {
            pNode = NULL;
            stack.pop();
            numPop += 1;
        }
    }

    m_parameters.pResults->updateShared(m_stats);
}   

void* GraphCompare::runThreaded(void* obj)
{
    GraphCompare compare(*reinterpret_cast<GraphCompareParameters*>(obj));
    compare.run();
    return NULL;
}

// Initialize the stack with a range of kmers
void GraphCompare::initializeStack(KmerRange range, GraphCompareStack& stack, const BWTVector& bwts, const BWTVector& rbwts)
{
    assert(range.start != -1 && range.end != -1);
    for(int64_t i = range.start; i < range.end; ++i)
    {
        // Decode i as a kmer value
        std::string str = range.decode(i);
        GraphCompareStackNode node;
        node.initialize(str, bwts, rbwts);
        stack.push(node);
    }
}

// Update the node and stack by either expanding to new nodes or updating
// the interval of the current node. Returns true if the node should
// be removed from the stack
bool GraphCompare::updateNodeAndStack(GraphCompareStackNode* pNode, GraphCompareStack& stack, const BWTVector& bwts, const BWTVector& rbwts)
{
    if(pNode->length == m_parameters.kmer)
        return true; // extend no further

    AlphaCount64 ext_count = pNode->getAggregateExtCount();
    bool hasExtension = ext_count.getNumNonZero() > 0;
    bool hasBranch = hasExtension && !ext_count.hasUniqueDNAChar();

    while(pNode->alphaIndex < DNA_ALPHABET::size)
    {
        char b = DNA_ALPHABET::getBase(pNode->alphaIndex);
        if(ext_count.get(b) == 0)
        {
            pNode->alphaIndex += 1;
            continue; // extension to b does not exist
        }

        if(hasBranch)
        {
            // Extension to b exists and there is a branch
            // Create a new element on the stack 
            GraphCompareStackNode branched = *pNode;
            branched.update(b, bwts, rbwts);
            pNode->alphaIndex += 1;
            stack.push(branched);

            // This node needs to be visited later for the other
            // possible extensions so it should not be popped
            return false;
        }
        else
        {
            // Update the interval for the current stack entry
            pNode->update(b, bwts, rbwts);
            return false;
        }
    }

    // No extensions were made, we are done with pNode
    return true;
}

// Returns true if the kmer represented by the node is a variant
bool GraphCompare::isVariantKmer(GraphCompareStackNode* pNode) const
{
    assert(pNode->length == m_parameters.kmer);
    return !pNode->intervalPairs[0].isValid() && pNode->intervalPairs[1].isValid();
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

    VariationBubbleBuilder builder;
    builder.setSourceIndex(bwts[varIndex], rbwts[varIndex]);
    builder.setTargetIndex(bwts[1 - varIndex], rbwts[1 - varIndex]);
    builder.setSourceString(str);

    //
    BubbleResult result = builder.run();
    if(result.returnCode == BRC_OK)
    {
        assert(!result.targetString.empty());
        assert(!result.sourceString.empty());
        markVariantSequenceKmers(result.sourceString);

        updateVariationCount(result);
        
        /*
        // Write to the variants file
        std::stringstream baseIDMaker;
        baseIDMaker << "base-" << 0;
        SeqItem item1 = { baseIDMaker.str(), result.targetString };
        item1.write(*m_pWriter);

        std::stringstream varIDMaker;
        varIDMaker << "variant-" << 0;
        SeqItem item2 = { varIDMaker.str(), result.sourceString };
        item2.write(*m_pWriter);
        */
    }

    // Update the results stats
    m_stats.numAttempted += 1;
    switch(result.returnCode)
    {
        case BRC_UNKNOWN:
            assert(false);
        case BRC_OK:
            m_stats.numBubbles += 1;
            break;
        case BRC_SOURCE_BROKEN:
            m_stats.numSourceBroken += 1;
            break;
        case BRC_SOURCE_BRANCH:
            m_stats.numSourceBranched += 1;
            break;
        case BRC_TARGET_BROKEN:
            m_stats.numTargetBroken += 1;
            break;
        case BRC_TARGET_BRANCH:
            m_stats.numTargetBranched += 1;
            break;
        case BRC_WALK_FAILED:
            m_stats.numWalkFailed += 1;
            break;
        case BRC_NO_SOLUTION:
            m_stats.numNoSolution += 1;
            break;
    }
    
    return result.returnCode == BRC_OK;
}

// Update the bit vector with the kmers that were assembled into str
void GraphCompare::markVariantSequenceKmers(const std::string& str)
{
    assert(str.size() >= m_parameters.kmer);
    size_t n = str.size() - m_parameters.kmer + 1;

    for(size_t i = 0; i < n; ++i)
    {
        std::string kseq = str.substr(i, m_parameters.kmer);
        BWTInterval interval = BWTAlgorithms::findInterval(m_parameters.pVariantBWT, kseq);
        if(interval.isValid())
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
                m_parameters.pBitVector->updateCAS(j, false, true);
        }

        // Mark the reverse complement k-mers too
        std::string rc_kseq = reverseComplement(kseq);
        interval = BWTAlgorithms::findInterval(m_parameters.pVariantBWT, rc_kseq);
        if(interval.isValid())
        {
            for(int64_t j = interval.lower; j <= interval.upper; ++j)
                m_parameters.pBitVector->updateCAS(j, false, true);
        }
    }
}

//
bool GraphCompare::isKmerMarked(const std::string& str) const
{
    assert(str.size() == m_parameters.kmer);
    BWTInterval interval = BWTAlgorithms::findInterval(m_parameters.pVariantBWT, str);
    if(interval.isValid())
    {
        if(m_parameters.pBitVector->test(interval.lower))
            return true;
    }
    
    // Mark the reverse complement k-mers too
    std::string rc_str = reverseComplement(str);
    interval = BWTAlgorithms::findInterval(m_parameters.pVariantBWT, rc_str);
    if(interval.isValid())
    {
        if(m_parameters.pBitVector->test(interval.lower))
            return true;
    }
    return false;
}

// Update the counts of each error type
void GraphCompare::updateVariationCount(const BubbleResult& result)
{
    std::string tsM;
    std::string vsM;
    StdAlnTools::makePaddedStrings(result.targetString, result.sourceString, tsM, vsM);

    assert(tsM.size() == vsM.size());
    
//    std::cout << "TSM: " << tsM << "\n";
//    std::cout << "VSM: " << vsM << "\n";

    // Count differences
    bool inIns = false;
    bool inDel = false;
    for(size_t i = 0; i < tsM.size(); ++i)
    {
        if(tsM[i] != '-' && vsM[i] != '-' && vsM[i] != tsM[i])
        {
            m_stats.numSubs += 1;
        }
        else if(tsM[i] == '-')
        {
            if(!inIns)
                m_stats.numInsertions += 1;
        }
        else if(vsM[i] == '-')
        {
            if(!inDel)
                m_stats.numDeletions += 1;
        }

        inIns = tsM[i] == '-';
        inDel = vsM[i] == '-';
    }
}

//
// GraphCompareAggregateResult
//
GraphCompareAggregateResults::GraphCompareAggregateResults()
{
    // Initialize mutex
    int ret = pthread_mutex_init(&m_mutex, NULL);
    if(ret != 0)
    {
        std::cerr << "Mutex initialization failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
GraphCompareAggregateResults::~GraphCompareAggregateResults()
{
    int ret = pthread_mutex_destroy(&m_mutex);
    if(ret != 0)
    {
        std::cerr << "Mutex destruction failed with error " << ret << ", aborting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

//
void GraphCompareAggregateResults::updateShared(const GraphCompareStats stats)
{
    pthread_mutex_lock(&m_mutex);
    m_stats.add(stats);
    pthread_mutex_unlock(&m_mutex);
}

//
void GraphCompareAggregateResults::printStats() const
{
    m_stats.print();
}
