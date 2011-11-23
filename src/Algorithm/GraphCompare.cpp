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
#include "LRAlignment.h"
#include "HapgenUtil.h"
#include "DindelRealignWindow.h"
#include "DindelUtil.h"
#include "HaplotypeBuilder.h"

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
    numHBFailed = 0;
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
    numHBFailed += other.numHBFailed;
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
    printf("Failed - haplotype builder failed: %d\n", numHBFailed);
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
    m_stats.clear();
}

//
GraphCompare::~GraphCompare()
{
}

GraphCompareResult GraphCompare::process(const SequenceWorkItem& item)
{
    GraphCompareResult result;
    SeqRecord currRead = item.read;
    std::string w = item.read.seq.toString();
    if(w.size() <= m_parameters.kmer)
    {
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
        BWTInterval interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pVariantBWT, 
                                                                    m_parameters.pVariantBWTCache, 
                                                                    kmer);

        // Check if this interval has been marked by a previous iteration of the loop
        assert(interval.isValid());
        if(m_parameters.pBitVector->test(interval.lower))
            continue;

        BWTInterval rc_interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pVariantBWT, 
                                                                       m_parameters.pVariantBWTCache, 
                                                                       reverseComplement(kmer));
        
        size_t count = interval.size();
        if(rc_interval.isValid())
            count += rc_interval.size();

        if(count >= m_parameters.kmerThreshold)
        {
            // Check if this k-mer is present in the other base index
            size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);

            if(base_count == 0)
            {
                // This is a variant kmer
                BWTVector bwts;
                bwts.push_back(m_parameters.pBaseBWT);
                bwts.push_back(m_parameters.pVariantBWT);
                //BubbleResult bubbleResult = processVariantKmer(kmer, count, bwts, 1);
                BubbleResult bubbleResult = processVariantKmerAggressive(kmer, count);

                if(bubbleResult.returnCode == BRC_OK)
                {
                    result.varStrings.push_back(bubbleResult.sourceString);
                    result.varCoverages.push_back(bubbleResult.sourceCoverage);

                    result.baseStrings.push_back(bubbleResult.targetString);
                    result.baseCoverages.push_back(bubbleResult.targetCoverage);
 
                    std::stringstream baseVCFSS;
                    std::stringstream variantVCFSS;
                    DindelReturnCode drc = DindelUtil::runDindelPair(bubbleResult.sourceString,
                                                                     bubbleResult.targetString,
                                                                     m_parameters,
                                                                     baseVCFSS,
                                                                     variantVCFSS);
                    
                    if(drc == DRC_OK)
                    {
                        std::cout << baseVCFSS.str() << "\n";
                        std::cout << variantVCFSS.str() << "\n";

                        result.baseVCFStrings.push_back(baseVCFSS.str());
                        result.variantVCFStrings.push_back(variantVCFSS.str());
                    }
                }
            }
        }

        // Update the bit vector
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
BubbleResult GraphCompare::processVariantKmer(const std::string& str, int count, const BWTVector& bwts, int varIndex)
{
    assert(varIndex == 0 || varIndex == 1);
    VariationBubbleBuilder builder;
    builder.setSourceIndex(bwts[varIndex]);
    builder.setTargetIndex(bwts[1 - varIndex]);
    builder.setSourceString(str, count);
    builder.setKmerThreshold(m_parameters.kmerThreshold);
    builder.setAllowedBranches(m_parameters.maxBranches);

    //
    BubbleResult result = builder.run();
    if(result.returnCode == BRC_OK)
    {
        assert(!result.targetString.empty());
        assert(!result.sourceString.empty());

        updateVariationCount(result);
        markVariantSequenceKmers(result.sourceString);
    }
    else
    {
        
        // Get all the kmers on this failed variant and mark them
        StringVector kmers = builder.getSourceKmers();
        for(size_t i = 0; i < kmers.size(); ++i)
            markVariantSequenceKmers(kmers[i]);
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
        case BRC_HB_FAILED:
            m_stats.numHBFailed += 1;
            break;
        case BRC_NO_SOLUTION:
            m_stats.numNoSolution += 1;
            break;
    }
    
    return result;
}

//
BubbleResult GraphCompare::processVariantKmerAggressive(const std::string& str, int /*count*/)
{

    std::cout << "VARIANT KMER: " << str << "\n";
    //
    std::string variant_str;
    size_t flanking_k_length = 0;
    bool found_variant_string = buildVariantStringConservative(str, variant_str, flanking_k_length);
    
    if(found_variant_string)
    {
        std::string transformed_str;
        bool found_transformed = transformVariantString(variant_str, transformed_str);
        if(found_transformed)
        {
            std::cout << "Transformed string!\n";
            StdAlnTools::globalAlignment(transformed_str, variant_str, true);
        }
    }

    BubbleResult result;
    result.returnCode = BRC_UNKNOWN;

    size_t hb_k = flanking_k_length;

    if(found_variant_string)
    {
        // Make a kmer count profile for the putative variant string

        std::cout << "VProfile(" << m_parameters.kmer << "):";
        IntVector countProfile = makeCountProfile(variant_str, m_parameters.kmer, m_parameters.pVariantBWT, 9);
        std::copy(countProfile.begin(), countProfile.end(), std::ostream_iterator<int>(std::cout, ""));
        std::cout << "\n";


        std::cout << "BProfile(" << m_parameters.kmer << "):";
        countProfile = makeCountProfile(variant_str, m_parameters.kmer, m_parameters.pBaseBWT, 9);
        std::copy(countProfile.begin(), countProfile.end(), std::ostream_iterator<int>(std::cout, ""));
        std::cout << "\n";

        std::cout << "BProfile(31): ";
        countProfile = makeCountProfile(variant_str, 31, m_parameters.pBaseBWT, 9);
        std::copy(countProfile.begin(), countProfile.end(), std::ostream_iterator<int>(std::cout, ""));
        std::cout << "\n";

        // Run haplotype builder
        // Set the start/end points of the haplotype builder
        std::string startAnchorSeq = variant_str.substr(0, hb_k); 
        std::string endAnchorSeq = variant_str.substr(variant_str.length() - hb_k);

        if(startAnchorSeq == endAnchorSeq)
            return result;

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
        builder.setIndex(m_parameters.pBaseBWT, NULL);
        builder.setKmerParameters(hb_k, m_parameters.bReferenceMode ? 1 : 2);

        // Run the builder
        HaplotypeBuilderReturnCode hbCode = builder.run();
        std::cout << "HBC: " << hbCode << "\n";
        HaplotypeBuilderResult hbResult;

        // The search was successful, build strings from the walks
        //BubbleResultCode bubbleCode = BRC_HB_FAILED;
        if(hbCode == HBRC_OK)
        {
            printf("Haplotype builder success\n");
            
            hbCode = builder.parseWalks(hbResult);
            if(hbCode == HBRC_OK)
            {
                std::string base_str = hbResult.haplotypes.front();
                StdAlnTools::globalAlignment(base_str, variant_str, true);
                result.returnCode = BRC_OK;
                result.targetString = base_str;
                result.sourceString = variant_str;
                markVariantSequenceKmers(variant_str);
            }
        }
    }

    return result;
}

//
bool GraphCompare::buildVariantStringGreedy(const std::string& startingKmer, std::string& outString, size_t& flanking_k_length)
{
    bool bLeftMatch = false;
    bool bRightMatch = false;

    size_t left_extensions = 0;
    size_t right_extensions = 0;

    // Iteratively extend the input kmer to the left and right until the ending kmer matches the reference
    size_t check_k = m_parameters.kmer;
    size_t extend_start_k = m_parameters.kmer;
    size_t extend_end_k = extend_start_k;
    size_t k_step = 5;
    size_t base_threshold = m_parameters.kmerThreshold;
    size_t extension_max = 200;

    assert(check_k <= m_parameters.kmer);

    outString = startingKmer;
    flanking_k_length = check_k;

    // Left extend
    while(!bLeftMatch)
    {
        if(left_extensions > extension_max)
            break;

        // Check if the first CHECK_K bases are found in the base index
        std::string checkMer = outString.substr(0, check_k);
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(checkMer, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);

        if(base_count >= base_threshold)
        {
            bLeftMatch = true;
            break;
        }
        else
        {
            // Extend the string by 1 base
            // Start with the highest kmer and relax the kmer until a match is found
            bool extended = false;
            for(size_t curr_k = extend_start_k; curr_k >= extend_end_k; curr_k -= k_step)
            {
                std::string kstr = outString.substr(0, curr_k);
                AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kstr, m_parameters.pVariantBWT, ED_ANTISENSE);
                if(extensionCounts.hasDNAChar())
                {
                    char b = extensionCounts.getMaxDNABase();
                    outString.insert(0, 1, b);
                    left_extensions += 1;
                    extended = true;
                    break;
                }
            }

            if(!extended)
                break;
        }
    }

    // Right extend
    while(!bRightMatch)
    {
        if(right_extensions > extension_max)
            break;

        // Check if the last check_l bases are found in the base index
        std::string checkMer = outString.substr(outString.length() - check_k, check_k);
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(checkMer, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);

        if(base_count >= base_threshold)
        {
            bRightMatch = true;
            break;
        }
        else
        {
            // Extend the string by 1 base
            // Start with the highest kmer and relax the kmer until a match is found
            bool extended = false;
            for(size_t curr_k = extend_start_k; curr_k >= extend_end_k; curr_k -= k_step)
            {
                std::string kstr = outString.substr(outString.length() - curr_k, curr_k);
                AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kstr, m_parameters.pVariantBWT, ED_SENSE);
                if(extensionCounts.hasDNAChar())
                {
                    char b = extensionCounts.getMaxDNABase();
                    outString.append(1, b);
                    right_extensions += 1;
                    extended = true;
                    break;
                }
            }

            if(!extended)
                break;
        }
    }

    if(bLeftMatch && bRightMatch)
        printf("Extension %s [%zu %zu] %s\n", bLeftMatch && bRightMatch ? "success" : "failed", left_extensions, right_extensions, outString.c_str());
    return bLeftMatch && bRightMatch;
}

//
bool GraphCompare::buildVariantStringConservative(const std::string& startingKmer, std::string& outString, size_t& flanking_k_length)
{
    bool bLeftMatch = false;
    bool bRightMatch = false;

    size_t left_extensions = 0;
    size_t right_extensions = 0;

    // Iteratively extend the input kmer to the left and right until the ending kmer matches the reference
    size_t check_k = m_parameters.kmer;
    size_t extend_start_k = m_parameters.kmer;
    size_t extend_end_k = 41;//extend_start_k;
    size_t k_step = 5;
    size_t base_threshold = m_parameters.bReferenceMode ? 1 : m_parameters.kmerThreshold;
    size_t extension_max = 200;

    assert(check_k <= m_parameters.kmer);

    outString = startingKmer;
    flanking_k_length = check_k;

    // Left extend
    while(!bLeftMatch)
    {
        if(left_extensions > extension_max)
            break;

        // Check if the first CHECK_K bases are found in the base index
        std::string checkMer = outString.substr(0, check_k);
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(checkMer, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);

        if(base_count >= base_threshold)
        {
            bLeftMatch = true;
            break;
        }

        // Extend the string by 1 base
        // Start with the highest kmer and relax the kmer until a match is found
        bool extended = false;
        for(size_t curr_k = extend_start_k; curr_k >= extend_end_k; curr_k -= k_step)
        {
            std::string kstr = outString.substr(0, curr_k);
            AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kstr, m_parameters.pVariantBWT, ED_ANTISENSE);

            // Search a unique extension
            char ext = '\0';
            int extCount = -1;
            bool bUniqueExt = extensionCounts.hasUniqueDNAChar();

            for(int i = 0; i < DNA_ALPHABET::size; ++i)
            {
                char b = DNA_ALPHABET::getBase(i);
                size_t ec = extensionCounts.get(b);
                if((ec > 0 && bUniqueExt) || ec >= m_parameters.kmerThreshold)
                {
                    extCount = ec;

                    // Check if we have set extIdx already
                    // If so this is a high-coverage branch and we stop
                    if(ext != '\0')
                    {
                        ext = '\0';
                        break;
                    }
                    ext = b;
                }
            }

            if(ext != '\0')
            {
                outString.insert(0, 1, ext);
                left_extensions += 1;
                extended = true;
                break;
            }
        }

        if(!extended)
            break;
    }

    // Right extend
    while(!bRightMatch)
    {
        if(right_extensions > extension_max)
            break;

        // Check if the last check_l bases are found in the base index
        std::string checkMer = outString.substr(outString.length() - check_k, check_k);
        size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(checkMer, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);

        if(base_count >= base_threshold)
        {
            bRightMatch = true;
            break;
        }

        // Extend the string by 1 base
        // Start with the highest kmer and relax the kmer until a match is found
        bool extended = false;
        for(size_t curr_k = extend_start_k; curr_k >= extend_end_k; curr_k -= k_step)
        {
            std::string kstr = outString.substr(outString.length() - curr_k, curr_k);
            AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kstr, m_parameters.pVariantBWT, ED_SENSE);

            // Search a unique extension
            char ext = '\0';
            int extCount = -1;
            bool bUniqueExt = extensionCounts.hasUniqueDNAChar();

            for(int i = 0; i < DNA_ALPHABET::size; ++i)
            {
                char b = DNA_ALPHABET::getBase(i);
                size_t ec = extensionCounts.get(b);
                if((ec > 0 && bUniqueExt) || ec >= m_parameters.kmerThreshold)
                {
                    extCount = ec;

                    // Check if we have set extIdx already
                    // If so this is a high-coverage branch and we stop
                    if(ext != '\0')
                    {
                        ext = '\0';
                        break;
                    }
                    ext = b;
                }
            }

            if(ext != '\0')
            {
                outString.append(1, ext);
                right_extensions += 1;
                extended = true;
                break;
            }
        }

        if(!extended)
            break;
    }

    if(bLeftMatch && bRightMatch)
        printf("Extension %s [%zu %zu] %s\n", bLeftMatch && bRightMatch ? "success" : "failed", left_extensions, right_extensions, outString.c_str());
    return bLeftMatch && bRightMatch;
}

// Transform inStr by substituting bases until all kmers covering it are found in the normal bwt
bool GraphCompare::transformVariantString(const std::string& inStr, std::string& outStr)
{
    size_t k = m_parameters.kmer;
    std::string inCopy = inStr;
    for(size_t i = k; i < inStr.size() - k + 1; ++i)
    {
        char cb = inStr[i];
        for(size_t j = 0; j < DNA_ALPHABET_SIZE; ++j)
        {
            char tb = DNA_ALPHABET::getBase(j);
            if(cb == tb)
                continue;
            inCopy[i] = tb;

            IntVector profile = makeCountProfile(inCopy, k, m_parameters.pBaseBWT, 100);
            bool allOK = true;
            for(size_t m = 0; m < profile.size(); ++m)
            {
                if(profile[m] < (int)m_parameters.kmerThreshold)
                {
                    allOK = false;
                    break;
                }
            }

            if(allOK)
            {
                outStr = inCopy;
                return true;
            }
        }

        inCopy[i] = cb;
    }

    return false;
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

IntVector GraphCompare::makeCountProfile(const std::string& str, size_t k, const BWT* pBWT, int max)
{
    IntVector out;
    if(str.size() < k)
        return out;

    for(size_t i = 0; i < str.size() - k + 1; ++i)
    {
        int count = BWTAlgorithms::countSequenceOccurrences(str.substr(i, k), pBWT);
        out.push_back(count > max ? max : count);
    }
    return out;
}

//
// GraphCompareAggregateResult
//
GraphCompareAggregateResults::GraphCompareAggregateResults(const std::string& fileprefix) : m_baseVCFFile(fileprefix + ".base.vcf","w"),
                                                                                            m_variantVCFFile(fileprefix + ".variant.vcf","w"),
                                                                                            m_numVariants(0)
{
    //
    m_pWriter = createWriter(fileprefix + ".strings.fa");
    
    //
    m_baseVCFFile.outputHeader("stub", "stub");
    m_variantVCFFile.outputHeader("stub", "stub");

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
    delete m_pWriter;

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

void GraphCompareAggregateResults::process(const SequenceWorkItem& /*item*/, const GraphCompareResult& result)
{
    assert(result.varStrings.size() == result.baseStrings.size());
    for(size_t i = 0; i < result.varStrings.size(); ++i)
    {
        // Write to the variants file
        std::stringstream baseIDMaker;
        std::stringstream baseMeta;
        baseIDMaker << "base-" << m_numVariants;
        baseMeta << "coverage=" << result.baseCoverages[i];

        SeqItem item1 = { baseIDMaker.str(), result.baseStrings[i] };
        item1.write(*m_pWriter, baseMeta.str());

        std::stringstream varIDMaker;
        std::stringstream varMeta;
        varIDMaker << "variant-" << m_numVariants;
        varMeta << "coverage=" << result.varCoverages[i];
        SeqItem item2 = { varIDMaker.str(), result.varStrings[i] };
        item2.write(*m_pWriter, varMeta.str());
        m_numVariants += 1;
    }

    assert(result.baseVCFStrings.size() == result.variantVCFStrings.size());
    for(size_t i = 0; i < result.baseVCFStrings.size(); ++i)
    {
        m_baseVCFFile.getOutputStream() << result.baseVCFStrings[i];
        m_variantVCFFile.getOutputStream() << result.variantVCFStrings[i];
    }
}

//
void GraphCompareAggregateResults::printStats() const
{
    m_stats.print();
}
