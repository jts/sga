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

GraphCompareResult GraphCompare::process(const SequenceWorkItem& item)
{
    GraphCompareResult result;
    SeqRecord currRead = item.read;
    std::string w = item.read.seq.toString();
    if(w.size() <= m_parameters.kmer)
    {
        result.temp = 1;
        return result;
    }

    // Perform a backwards search using the read sequence
    // Check which k-mers have already been visited using the
    // shared bitvector. If any bit in the range [l,u] is set
    // for a suffix of the read, then we do not visit those kmers
    // later
    int len = w.size();
    int num_kmers = len - m_parameters.kmer + 1;
    std::vector<bool> visitedKmers(false, num_kmers);

    int j = len - 1;
    char curr = w[j];
    BWTInterval interval;
    BWTAlgorithms::initInterval(interval, curr, m_parameters.pVariantBWT);
    --j;

    for(;j >= 0; --j)
    {
        curr = w[j];
        BWTAlgorithms::updateInterval(interval, curr, m_parameters.pVariantBWT);
        assert(interval.isValid());

        // At this point interval represents the suffix [j,len)
        // Check if the starting point of this interval is set
        if(j < num_kmers)
            visitedKmers[j] = m_parameters.pBitVector->test(interval.lower);
    }
    
    // Process the kmers that have not been previously visited
    for(j = 0; j < num_kmers; ++j)
    {
        if(visitedKmers[j])
            continue; // skip
        std::string kmer = w.substr(j, m_parameters.kmer);
        
        // Get the interval for this kmer
        BWTInterval interval = BWTAlgorithms::findInterval(m_parameters.pVariantBWT, kmer);
        BWTInterval rc_interval = BWTAlgorithms::findInterval(m_parameters.pVariantBWT, reverseComplement(kmer));

        // Mark the bit vector so it wont be visited again
        assert(interval.isValid());
        
        size_t count = interval.size();
        if(rc_interval.isValid())
            count += rc_interval.size();

        if(count > m_parameters.kmerThreshold)
        {
            // Check if this k-mer is present in the other base index
            size_t base_count = BWTAlgorithms::countSequenceOccurrences(kmer, m_parameters.pBaseBWT);
            if(base_count == 0)
            {
                // This is a variant kmer
                BWTVector bwts;
                bwts.push_back(m_parameters.pBaseBWT);
                bwts.push_back(m_parameters.pVariantBWT);

                BWTVector rbwts;
                rbwts.push_back(m_parameters.pBaseRevBWT);
                rbwts.push_back(m_parameters.pVariantRevBWT);
                processVariantKmer(kmer, bwts, rbwts, 1);
            }
        }

        for(int64_t i = interval.lower; i <= interval.upper; ++i)
            m_parameters.pBitVector->set(i, true);

        for(int64_t i = rc_interval.lower; i <= rc_interval.upper; ++i)
            m_parameters.pBitVector->set(i, true);
        
    }
        
    return result;
}

void GraphCompare::updateSharedStats(GraphCompareAggregateResults* pSharedStats)
{
    pSharedStats->updateShared(m_stats);
    m_stats.clear();
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

void GraphCompareAggregateResults::process(const SequenceWorkItem& item, const GraphCompareResult& result)
{
    (void)item;
    (void)result;
}

//
void GraphCompareAggregateResults::printStats() const
{
    m_stats.print();
}
