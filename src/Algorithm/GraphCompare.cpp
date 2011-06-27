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

    //str.append(1,b);
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

    //str.append(1,b);
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


}

//
GraphCompare::~GraphCompare()
{

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
    size_t numPush = 0;
    size_t numPop = 0;

    while(!stack.empty())
    {

        GraphCompareStackNode* pNode = &stack.top();
        
        // Progress counter bookkeeping
        loops += 1;
        if(stack.size() > maxStack)
            maxStack = stack.size();

        if(loops % 1000000 == 0)
        {
            std::cout << "Loop: " << loops << "\n";
            std::cout << "Kmers: " << kmersFound << "\n";
            pNode->print();
        }

        bool doPop = true;

        if(pNode->length == (int)m_kmer)
        {
            // do something
            kmersFound += 1;

            if(pNode->intervalPairs[1].isValid() && !pNode->intervalPairs[0].isValid())
                variantKmersFound += 1;
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

}   
