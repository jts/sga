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
    m_pWriter = createWriter("variants.fa");

    m_numBubbles = 0;
    m_numAttempted = 0;
    m_numTargetBranched = 0;
    m_numSourceBranched = 0;
    m_numTargetBroken = 0;
    m_numSourceBroken = 0;
    m_numWalkFailed = 0;
    m_numNoSolution = 0;
        
    m_numInsertions = 0;
    m_numDeletions = 0;
    m_numSubs = 0;
}

//
GraphCompare::~GraphCompare()
{
    delete m_pUsedVariantKmers;
    delete m_pWriter;
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

    printf("Total attempts: %d\n", m_numAttempted);
    printf("Total bubbles: %d\n", m_numBubbles);
    printf("Failed - target branched: %d\n", m_numTargetBranched);
    printf("Failed - source branched: %d\n", m_numSourceBranched);
    printf("Failed - target broken: %d\n", m_numTargetBroken);
    printf("Failed - source broken: %d\n", m_numSourceBroken);
    printf("Failed - no walk: %d\n", m_numWalkFailed);
    printf("Failed - no solution: %d\n", m_numNoSolution);
    
    printf("Num subs found: %d\n", m_numSubs);
    printf("Num insertions found: %d\n", m_numInsertions);
    printf("Num deletions found: %d\n", m_numDeletions);
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

        // Write to the variants file
        std::stringstream baseIDMaker;
        baseIDMaker << "base-" << m_numBubbles;
        SeqItem item1 = { baseIDMaker.str(), result.targetString };
        item1.write(*m_pWriter);

        std::stringstream varIDMaker;
        varIDMaker << "variant-" << m_numBubbles;
        SeqItem item2 = { varIDMaker.str(), result.sourceString };
        item2.write(*m_pWriter);
    }

    // Update the results stats
    m_numAttempted += 1;
    switch(result.returnCode)
    {
        case BRC_UNKNOWN:
            assert(false);
        case BRC_OK:
            m_numBubbles += 1;
            break;
        case BRC_SOURCE_BROKEN:
            m_numSourceBroken += 1;
            break;
        case BRC_SOURCE_BRANCH:
            m_numSourceBranched += 1;
            break;
        case BRC_TARGET_BROKEN:
            m_numTargetBroken += 1;
            break;
        case BRC_TARGET_BRANCH:
            m_numTargetBranched += 1;
            break;
        case BRC_WALK_FAILED:
            m_numWalkFailed += 1;
            break;
        case BRC_NO_SOLUTION:
            m_numNoSolution += 1;
            break;
    }
    
    return result.returnCode == BRC_OK;
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
            m_numSubs += 1;
        }
        else if(tsM[i] == '-')
        {
            if(!inIns)
            {
                m_numInsertions += 1;
            }
        }
        else if(vsM[i] == '-')
        {
            if(!inDel)
            {
                m_numDeletions += 1;
            }
        }

        if(tsM[i] == '-')
            inIns = true;
        else
            inIns = false;
        if(vsM[i] == '-')
            inDel = true;
        else
            inDel = false;
    }
}

