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
#include "ReadCoherentHaplotypeBuilder.h"
#include "OverlapHaplotypeBuilder.h"
#include "BuilderCommon.h"
#include "Profiler.h"


// #define GRAPH_DIFF_DEBUG 1

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
    PROFILE_FUNC("GraphCompare::process")
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
    // We only try to find variants from one k-mer per read. 
    // This heurestic handles the situation where we process a kmer,
    // spend a lot of time trying to assemble it into a string and fail,
    // then try again with the next k-mer in the read. Since all the k-mers
    // in a read are connected in the de Bruijn graph, if the first attempt
    // to assemble a k-mer into a variant fails, it is very unlikely that
    // subsequent attempts will succeed.
    bool variantAttempted = false;
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

        bool both_directions = interval.size() > 0 && rc_interval.size() > 0;

        if(count >= m_parameters.kmerThreshold && count < m_parameters.maxKmerThreshold && !variantAttempted && both_directions)
        {
            // Check if this k-mer is present in the other base index
            size_t base_count = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);

            if(base_count == 0)
            {
                // This is a variant kmer
                variantAttempted = true;
                BWTVector bwts;
                bwts.push_back(m_parameters.pBaseBWT);
                bwts.push_back(m_parameters.pVariantBWT);
                std::cout << "Variant read: " << w << "\n";
                //BubbleResult bubbleResult = processVariantKmer(kmer, count, bwts, 1);
                GraphBuildResult build_result = processVariantKmerAggressive(kmer, count);

                // Mark the kmers of the variant haplotypes as being visited
                for(size_t vhi = 0; vhi < build_result.variant_haplotypes.size(); ++vhi)
                    markVariantSequenceKmers(build_result.variant_haplotypes[vhi]);
                
                //printf("Build %zu var haps %zu base haps\n", variant_haplotypes.size(), base_haplotypes.size());

                if(build_result.variant_haplotypes.size() > 0 /*&& build_result.base_haplotypes.size() > 0*/)
                {
                    std::cout << "Running dindel\n";
                    std::stringstream baseVCFSS;
                    std::stringstream variantVCFSS;
                    
                    DindelReturnCode drc = DindelUtil::runDindelPairMatePair(kmer,
                                                                             build_result.base_haplotypes,
                                                                             build_result.variant_haplotypes,
                                                                             m_parameters,
                                                                             baseVCFSS,
                                                                             variantVCFSS);
                    
                    std::cout << "Dindel returned " << drc << "\n";
                    std::cout << "base: " << baseVCFSS.str() << "\n";
                    std::cout << "vari: " << variantVCFSS.str() << "\n";
			    
                    if(drc == DRC_OK)
                    {                        
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
GraphBuildResult GraphCompare::processVariantKmerAggressive(const std::string& str, int count)
{
    PROFILE_FUNC("GraphCompare::processVariantKmerAggressive")

#ifdef GRAPH_DIFF_DEBUG
    std::cout << "Processing variant kmer " << str << " with depth: " << count << "\n";
#endif
    (void)count;

    size_t haplotype_builder_kmer = m_parameters.kmer;

    //
    GraphBuildResult result;

    /*
    ReadCoherentHaplotypeBuilder rc_builder;
    rc_builder.setInitialHaplotype(str);
    rc_builder.setParameters(m_parameters);
    rc_builder.run(result.variant_haplotypes);
    */

    OverlapHaplotypeBuilder overlap_builder(m_parameters);
    overlap_builder.setInitialHaplotype(str);
    overlap_builder.run(result.variant_haplotypes);

    // Haplotype QC
    // Calculate the maximum k such that every kmer is present in the variant and base BWT
    // The difference between these values must be at least MIN_COVER_K_DIFF
    size_t MIN_COVER_K_DIFF = 20;
    StringVector temp_haplotypes;
    for(size_t i = 0; i < result.variant_haplotypes.size(); ++i)
    {
        size_t max_variant_k = calculateMaxCoveringK(result.variant_haplotypes[i], 1, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache);
        size_t max_base_k = calculateMaxCoveringK(result.variant_haplotypes[i], 1, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache);
        printf("MVK: %zu MBK: %zu\n", max_variant_k, max_base_k);
        if( max_variant_k > max_base_k && max_variant_k - max_base_k >= MIN_COVER_K_DIFF && max_base_k < 21)
            temp_haplotypes.push_back(result.variant_haplotypes[i]);
    }
    result.variant_haplotypes.swap(temp_haplotypes);

    bool found_variant_string = result.variant_haplotypes.size() > 0;

    for(size_t i = 0; i < result.variant_haplotypes.size(); ++i)
        printf("Assembly[%zu]: %s\n", i, result.variant_haplotypes[i].c_str());

    // JTS: disable base string generation which can lead to false calls if wrong
    if(found_variant_string &&  false)
    {
        // Run haplotype builder on the normal graph
        for(size_t i = 0; i < result.variant_haplotypes.size(); ++i)
        {
            const std::string& current_variant_haplotype = result.variant_haplotypes[i];

            std::string startAnchorSeq = current_variant_haplotype.substr(0, haplotype_builder_kmer); 
            std::string endAnchorSeq = current_variant_haplotype.substr(current_variant_haplotype.length() - haplotype_builder_kmer);

#ifdef GRAPH_DIFF_DEBUG
            // Make a kmer count profile for the putative variant string
            std::cout << "VProfile(" << m_parameters.kmer << "):";
            IntVector countProfile = makeCountProfile(current_variant_haplotype, m_parameters.kmer, m_parameters.pVariantBWT, 9);
            std::copy(countProfile.begin(), countProfile.end(), std::ostream_iterator<int>(std::cout, ""));
            std::cout << "\n";


            std::cout << "BProfile(" << m_parameters.kmer << "):";
            countProfile = makeCountProfile(current_variant_haplotype, m_parameters.kmer, m_parameters.pBaseBWT, 9);
            std::copy(countProfile.begin(), countProfile.end(), std::ostream_iterator<int>(std::cout, ""));
            std::cout << "\n";

            std::cout << "BProfile(31): ";
            countProfile = makeCountProfile(current_variant_haplotype, 31, m_parameters.pBaseBWT, 9);
            std::copy(countProfile.begin(), countProfile.end(), std::ostream_iterator<int>(std::cout, ""));
            std::cout << "\n";


            std::cout << "Mapping locations for start: " << startAnchorSeq << "\n";
            showMappingLocations(startAnchorSeq);
            std::cout << "Mapping locations for end: " << endAnchorSeq << "\n";
            showMappingLocations(endAnchorSeq);
#endif

            if(startAnchorSeq == endAnchorSeq)
                continue; // degenerate sequence, skip

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
            builder.setKmerParameters(haplotype_builder_kmer, m_parameters.bReferenceMode ? 1 : 2);

            // Run the builder
            HaplotypeBuilderReturnCode hbCode = builder.run();
//#ifdef GRAPH_DIFF_DEBUG
            std::cout << "HBC: " << hbCode << " VS: " << current_variant_haplotype << "\n";
//#endif
            HaplotypeBuilderResult hbResult;

            // The search was successful, build strings from the walks
            if(hbCode == HBRC_OK) {
    
                hbCode = builder.parseWalks(hbResult);
                if(hbCode == HBRC_OK) {
                    result.base_haplotypes.insert(result.base_haplotypes.end(), hbResult.haplotypes.begin(), hbResult.haplotypes.end());
                }
            }
        }
    }
#ifdef GRAPH_DIFF_DEBUG
    else
    {
        std::cout << "Variant string not assembled\n";
    }
#endif

    return result;
}

// Graph-based method of building a variant string until it meets the graph of the normal
bool GraphCompare::buildVariantStringGraph(const std::string& startingKmer, StringVector& haplotypes)
{
    PROFILE_FUNC("GraphCompare::buildVariantStringGraph")

    std::map<std::string, int> kmerCountMap;

    // We search until we find the first common vertex in each direction

    size_t MIN_TARGET_COUNT = m_parameters.bReferenceMode ? 1 : 2;
    size_t MAX_ITERATIONS = 2000;
    size_t MAX_SIMULTANEOUS_BRANCHES = 20;
    size_t MAX_TOTAL_BRANCHES = 50;

    // Tracking stats
    size_t max_simul_branches_used = 0;
    size_t total_branches = 0;
    size_t iterations = 0;

    // Initialize the graph
    StringGraph* pGraph = new StringGraph;
    BuilderExtensionQueue queue;

    Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(startingKmer, startingKmer);
    pVertex->setColor(GC_BLACK);
    pGraph->addVertex(pVertex);

    // Add the vertex to the extension queue
    queue.push(BuilderExtensionNode(pVertex, ED_SENSE));
    queue.push(BuilderExtensionNode(pVertex, ED_ANTISENSE));

    std::vector<Vertex*> sense_join_vector;
    std::vector<Vertex*> antisense_join_vector;

    // Perform the extension. The while conditions are heuristics to avoid searching
    // the graph too much 
    while(!queue.empty() && iterations++ < MAX_ITERATIONS && queue.size() < MAX_SIMULTANEOUS_BRANCHES && total_branches < MAX_TOTAL_BRANCHES)
    {
        if(queue.size() > max_simul_branches_used)
            max_simul_branches_used = queue.size();

        BuilderExtensionNode curr = queue.front();
        queue.pop();

        /*
        // We have found a join in this direction, do not continue
        if(joinFound[curr.direction])
            continue;
        */

        // Calculate de Bruijn extensions for this node
        std::string vertStr = curr.pVertex->getSeq().toString();
        AlphaCount64 extensionCounts = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(vertStr, m_parameters.pVariantBWT, curr.direction);

        std::string extensionsUsed;
        for(size_t i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            size_t count = extensionCounts.get(b);
            //bool acceptExt = count >= m_parameters.kmerThreshold || (count > 0 && extensionCounts.hasUniqueDNAChar());
            bool acceptExt = count >= m_parameters.minKmerThreshold;
            if(!acceptExt)
                continue;

            extensionsUsed.push_back(b);
            std::string newStr = BuilderCommon::makeDeBruijnVertex(vertStr, b, curr.direction);
            kmerCountMap[newStr] = count;

            // Create the new vertex and edge in the graph
            // Skip if the vertex already exists
            if(pGraph->getVertex(newStr) != NULL)
                continue;
            
            // Allocate the new vertex and add it to the graph
            Vertex* pVertex = new(pGraph->getVertexAllocator()) Vertex(newStr, newStr);
            pVertex->setColor(GC_BLACK);
            pGraph->addVertex(pVertex);

            // Add edges
            BuilderCommon::addSameStrandDeBruijnEdges(pGraph, curr.pVertex, pVertex, curr.direction);
            
            // Check if this sequence is present in the FM-index of the target
            // If so, it is the join point of the de Bruijn graph and we extend no further.
            size_t targetCount = BWTAlgorithms::countSequenceOccurrences(newStr, m_parameters.pBaseBWT);
            if(targetCount >= MIN_TARGET_COUNT)
            {
                if(curr.direction == ED_SENSE)
                    sense_join_vector.push_back(pVertex);
                else
                    antisense_join_vector.push_back(pVertex);
            }
            else
            {
                // Add the vertex to the extension queue
                queue.push(BuilderExtensionNode(pVertex, curr.direction));
            }
        }
        
        // Update the total number of times we branches the search
        if(!extensionsUsed.empty())
            total_branches += extensionsUsed.size() - 1;
    }

    /*
    std::string result_str = (pSenseJoin != NULL && pAntisenseJoin != NULL) ? "OK" : "FAIL";
    printf("VariantStringGraph\t%s\tMS:%zu\tTB:%zu\tNI:%zu\n", result_str.c_str(), max_simul_branches_used, total_branches, iterations);
    */

    // If the graph construction was successful, walk the graph
    // between the endpoints to make a string
    // Generate haplotypes between every pair of antisense/sense join vertices
    for(size_t i = 0; i < antisense_join_vector.size(); ++i) {
        for(size_t j = 0; j < sense_join_vector.size(); ++j) {
            SGWalkVector outWalks;
            SGSearch::findWalks(antisense_join_vector[i],
                                sense_join_vector[j],
                                ED_SENSE,
                                100000, // max distance to search
                                10000, // max nodes to search
                                true, // exhaustive search
                                outWalks);

            for(size_t k = 0; k < outWalks.size(); ++k)
                haplotypes.push_back(outWalks[k].getString(SGWT_START_TO_END));
        }
    }
    
    delete pGraph;
    return !haplotypes.empty();
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
size_t GraphCompare::calculateMaxCoveringK(const std::string& sequence, int min_depth, const BWT* pBWT, const BWTIntervalCache* pBWTCache)
{
    size_t k = 99;
    size_t min_k = 15;
    while(k >= min_k)
    {
        if(sequence.size() < k)
            continue;
        bool covered = true;
        size_t nk = sequence.size() - k + 1;
        for(size_t i = 0; i < nk; ++i)
        {
            std::string kmer = sequence.substr(i, k);
            int c = BWTAlgorithms::countSequenceOccurrencesWithCache(kmer, pBWT, pBWTCache);

            if(c < min_depth)
            {
                covered = false;
                break;
            }
        }

        if(covered)
            return k;
        k -= 1;
    }

    return 0;
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

void GraphCompare::debug(const std::string& debugFilename)
{
        (void)debugFilename;
#if 0
    std::cout << "Debug file: " << debugFilename << "\n";   
    std::istream* pReader = createReader(debugFilename);
    
    std::string line;
    StringVector kmerVector;
    while(getline(*pReader, line))
    {
        StringVector fields = split(line, '\t');

        if(fields[0] == "NEXT")
        {
            // Process the set of kmers
            for(size_t i = 0; i < kmerVector.size(); ++i)
            {
                std::string kmer = kmerVector[i];
                size_t bc = BWTAlgorithms::countSequenceOccurrences(kmer, m_parameters.pBaseBWT);
                size_t vc = BWTAlgorithms::countSequenceOccurrences(kmer, m_parameters.pVariantBWT);
                size_t rc = BWTAlgorithms::countSequenceOccurrences(kmer, m_parameters.pReferenceBWT);
                AlphaCount64 aec = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kmer, m_parameters.pVariantBWT, ED_ANTISENSE);
                AlphaCount64 sec = BWTAlgorithms::calculateDeBruijnExtensionsSingleIndex(kmer, m_parameters.pVariantBWT, ED_SENSE);
                std::cout << kmer << " BC:" << bc << " VC:" << vc << " RC:" << rc << " SEC:" << sec << " AEC: " << aec << "\n";
            }
         
            if(!kmerVector.empty())
            {
                BubbleResult bubbleResult = processVariantKmerAggressive(kmerVector[kmerVector.size() / 2], 0);

                if(bubbleResult.returnCode == BRC_OK)
                {
                    std::stringstream baseVCFSS;
                    std::stringstream variantVCFSS;
                    std::string kmer = ".";
                    DindelReturnCode drc = DindelUtil::runDindelPairMatePair(kmer,
                                                                             bubbleResult.sourceString,
                                                                             bubbleResult.targetString,
                                                                             m_parameters,
                                                                             baseVCFSS,
                                                                             variantVCFSS);
                    
                    if(drc == DRC_OK)
                    {
                        std::cout << "BASE: " << baseVCFSS.str() << "\n";
                        std::cout << "VAR: " << variantVCFSS.str() << "\n";
                    }
                    else
                    {
                        std::cout << "Dindel fail: " << drc << "\n";
                    }
                }
            }
            kmerVector.clear();
        }
        else if(fields.size() == 4)
        {
            kmerVector.push_back(fields[0]);
        }
        else
        {
            std::cout << "Processing " << line << "\n";
        }
    }
#endif
    
}

// test kmers from a file
void GraphCompare::testKmersFromFile(const std::string& kmerFilename)
{
    std::cout << "Kmer file: " << kmerFilename << "\n";
    std::istream* pReader = createReader(kmerFilename);

    std::string line;
    while(getline(*pReader, line))
    {
        std::cout << "Processing " << line << "\n";
        StringVector fields = split(line, '\t');
        std::string kmer = fields[0];

        if(kmer.size() != m_parameters.kmer)
        {
            std::cout << "Incorrect kmer size: " << kmer << " kmer size: " << kmer.size() << " kmer length parameter: " << m_parameters.kmer << "\n";
            continue;
        }
        testKmer(kmer);
    }
}

void GraphCompare::testKmer(const std::string& kmer)
{
    int count = 1;
    GraphBuildResult build_result = processVariantKmerAggressive(kmer, count);

    if(build_result.variant_haplotypes.size() > 0 && build_result.base_haplotypes.size() > 0)
    {
        std::cout << "Haplotypes successfully built. Aligning first pair.\n";
        StdAlnTools::globalAlignment(build_result.base_haplotypes.front(), build_result.variant_haplotypes.front(), true);

        std::stringstream baseVCFSS;
        std::stringstream variantVCFSS;

        DindelReturnCode drc = DindelUtil::runDindelPairMatePair(kmer,
                                                                 build_result.base_haplotypes,
                                                                 build_result.variant_haplotypes,
                                                                 m_parameters,
                                                                 baseVCFSS,
                                                                 variantVCFSS);
        
        std::cout << "base:    " << baseVCFSS.str() << "\n";
        std::cout << "variant: " << variantVCFSS.str() << "\n";
        
        if(drc == DRC_OK)
            std::cout << "DINDEL says: OK.\n";
        else
            std::cout << "DINDEL error: " << drc <<"\n";
    } 
    else 
    {
        
        std::cout << "Error. Not enough haplotypes VH:" << build_result.variant_haplotypes.size() << " BH: " << build_result.base_haplotypes.size() << "\n";
    }
}

void GraphCompare::showMappingLocations(const std::string& str)
{
    BWTInterval ref_interval = BWTAlgorithms::findInterval(m_parameters.pReferenceBWT, str);
    for(int64_t i = ref_interval.lower; i <= ref_interval.upper; ++i)
        std::cout << "Location: " << m_parameters.pReferenceSSA->calcSA(i, m_parameters.pReferenceBWT) << "\n";
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
