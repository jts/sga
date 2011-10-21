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
GraphCompare::GraphCompare(const GraphCompareParameters& params) : m_parameters(params), 
                                                                   m_baseVCFFile("base.vcf","w"),
                                                                   m_variantVCFFile("variant.vcf","w")
{
    m_baseVCFFile.outputHeader("stub", "stub");
    m_variantVCFFile.outputHeader("stub", "stub");
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
        BWTInterval interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, kmer);

        // Check if this interval has been marked by a previous iteration of the loop
        assert(interval.isValid());
        if(m_parameters.pBitVector->test(interval.lower))
            continue;

        BWTInterval rc_interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pVariantBWT, m_parameters.pVariantBWTCache, reverseComplement(kmer));
        
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

                BubbleResult bubbleResult = processVariantKmer(kmer, count, bwts, 1);

                if(bubbleResult.returnCode == BRC_OK)
                {
                    result.varStrings.push_back(bubbleResult.sourceString);
                    result.varCoverages.push_back(bubbleResult.sourceCoverage);

                    result.baseStrings.push_back(bubbleResult.targetString);
                    result.baseCoverages.push_back(bubbleResult.targetCoverage);
 
                    runDindelIndividual(bubbleResult.sourceString, bubbleResult.targetString);
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

// Perform the Dindel realignment and inference framework
void GraphCompare::runDindelIndividual(const std::string& normalString, const std::string& variantString)
{
    StringVector inHaplotypes;
    inHaplotypes.push_back(normalString);
    inHaplotypes.push_back(variantString);

    // Extract the reads from the normal and variant data sets that match each haplotype
    
    // Normal reads
    SeqItemVector normalReads;
    SeqItemVector normalReadMates;
    SeqItemVector normalRCReads;
    SeqItemVector normalRCReadMates;
        
    HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache,
                                      m_parameters.pBaseSSA, m_parameters.kmer, false, &normalReads, &normalReadMates);

    HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache,
                                      m_parameters.pBaseSSA, m_parameters.kmer, true, &normalRCReads, &normalRCReadMates);

    // Variant reads
    SeqItemVector variantReads;
    SeqItemVector variantReadMates;
    SeqItemVector variantRCReads;
    SeqItemVector variantRCReadMates;

    HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache,
                                      m_parameters.pVariantSSA, m_parameters.kmer, false, &variantReads, &variantReadMates);

    HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache,
                                      m_parameters.pVariantSSA, m_parameters.kmer, true, &variantRCReads, &variantRCReadMates);

    // Align the haplotypes to the reference genome to generate candidate alignments
    HapgenAlignmentVector candidateAlignments;
    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        HapgenUtil::alignHaplotypeToReference(inHaplotypes[i], m_parameters.pReferenceBWT, 
                                              m_parameters.pReferenceSSA, candidateAlignments);
    }

    // Remove duplicate or bad alignment pairs
    HapgenUtil::coalesceAlignments(candidateAlignments);

    // Score each candidate alignment against the mates of all the variant reads
    int bestCandidate = -1;
    double bestAverageScoreFrac = 0.0f;
    double secondBest = 0.0f;
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        // Compute the average score of the reads' mates to the flanking sequence
        StringVector referenceFlanking;
        std::cout << "Processing " << candidateAlignments[i] << "\n";
        HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i], m_parameters.pRefTable, 
                                           1000, inHaplotypes, referenceFlanking);


        // If valid flanking haplotypes could not be made, skip this alignment
        if(referenceFlanking.empty())
            continue;

        // Realign the mates
        LocalAlignmentResultVector localAlignments = HapgenUtil::alignReadsLocally(referenceFlanking[0], variantReadMates);
        LocalAlignmentResultVector localAlignmentsRC = HapgenUtil::alignReadsLocally(referenceFlanking[0], variantRCReadMates);

        // Merge alignments
        localAlignments.insert(localAlignments.end(), localAlignmentsRC.begin(), localAlignmentsRC.end());

        double sum = 0.0f;
        double count = 0.0f;

        for(size_t j = 0; j < localAlignments.size(); ++j)
        {
            double max_score = localAlignments[j].queryEndPosition - localAlignments[j].queryStartPosition;
            double frac = (double)localAlignments[j].score / max_score;
            //printf("Score: %d frac: %lf\n", localAlignments[j].score, frac);
            sum += frac;
            count += 1;
        }

        double score = sum / count;
        if(score > bestAverageScoreFrac)
        {
            secondBest = bestAverageScoreFrac;
            bestAverageScoreFrac = score;
            bestCandidate = i;
        }
        else if(score > secondBest)
        {
            secondBest = score;
        }

        printf("Alignment %zu mate-score: %lf\n", i, score);
    }

    printf("total alignments: %zu best score: %lf\n", candidateAlignments.size(), bestAverageScoreFrac);

    if(bestCandidate == -1)
    {
        std::cout << "No good alignment candidate\n";
        return;
    }

    if(bestAverageScoreFrac < 0.9f)
    {
        std::cout << "Skipping marginal pair alignment\n";
        return;
    }

    if(bestAverageScoreFrac - secondBest < 0.05f)
    {
        std::cout << "Best score too close to second best\n";
        return;
    }

    //
    // Run dindel against the haplotypes for the best alignment position
    // We do this for the normal and variant reads separately
    //

    // Generate the input haplotypes for dindel
    int FLANKING_SIZE = 0;
    StringVector flankingHaplotypes;
    HapgenUtil::makeFlankingHaplotypes(candidateAlignments[bestCandidate], m_parameters.pRefTable, 
                                       FLANKING_SIZE, inHaplotypes, flankingHaplotypes);

    // We need at least 2 haplotypes
    assert(flankingHaplotypes.size() >= 2);
    assert(flankingHaplotypes[0].size() > 0);

    //
    // Run Dindel
    //
    double MAP_QUAL = 40.0;
    int BASE_QUAL = 20;

    for(size_t i = 0; i <= 1; ++i)
    {
        SeqItemVector& fwdReads = (i == 0) ? normalReads : variantReads;
        SeqItemVector& rcReads = (i == 0) ? normalRCReads : variantRCReads;

        // Create dindel reads for the normal reads
        std::vector<DindelRead> dReads;
        for(size_t j = 0; j < fwdReads.size(); ++j) 
            dReads.push_back(DindelRead(fwdReads[j], std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, true));

        for(size_t j = 0; j < rcReads.size(); ++j)
        {
            rcReads[j].seq.reverseComplement();
            dReads.push_back(DindelRead(rcReads[j], std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, false));
        }

        int dindelRefStart = candidateAlignments[bestCandidate].position + 1; // VCF coordinates are 1-based
        std::stringstream refName;
        refName << m_parameters.pRefTable->getRead(candidateAlignments[bestCandidate].referenceID).id;

        std::string dindelRef = flankingHaplotypes[0]; // First flanking haplotype is of the reference
        StringVector nonReference(flankingHaplotypes.begin()+1, flankingHaplotypes.end());
        std::cout << "Running dindel on " << nonReference.size() << " haplotypes and " << dReads.size() << " reads\n";

        try
        {
            DindelWindow dWindow(nonReference, dindelRef, dindelRefStart, refName.str() );

            DindelRealignParameters dRealignParameters("addSNPMaxSNPs:0");
            DindelRealignWindow dRealignWindow(&dWindow, dReads, dRealignParameters);

            dRealignWindow.run("hmm", i == 0 ? m_baseVCFFile : m_variantVCFFile);
        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            exit(EXIT_FAILURE);
        }
    }
}

// Perform the Dindel realignment and inference framework
void GraphCompare::runDindelFull(const std::string& normalString, const std::string& variantString)
{
    StringVector inHaplotypes;
    inHaplotypes.push_back(normalString);
    inHaplotypes.push_back(variantString);

    // First, calculate alignments of the input haplotypes to the reference
    HapgenAlignmentVector candidateAlignments;
    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        HapgenUtil::alignHaplotypeToReference(inHaplotypes[i],
                                              m_parameters.pReferenceBWT, 
                                              m_parameters.pReferenceSSA,
                                              candidateAlignments);
    }


    // Remove duplicate or bad alignment pairs
    HapgenUtil::coalesceAlignments(candidateAlignments);
    std::cout << "Found " << candidateAlignments.size() << " alignments for all strings\n";

    // Join each haplotype with flanking sequence from the reference genome for each alignment
    // This function also adds a haplotype (with flanking sequence) for the piece of the reference
    bool success = true;
    int FLANKING_SIZE = 500;
    StringVector flankingHaplotypes;
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        std::cout << "Score: " << candidateAlignments[i].score << "\n";
        success = HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i], 
                                                     m_parameters.pRefTable, 
                                                     FLANKING_SIZE, 
                                                     inHaplotypes,
                                                     flankingHaplotypes);
        if(!success)
            break;
    }

    // Extract all read pairs that match a k-mer to any haplotype
    if(candidateAlignments.size() > 0 && flankingHaplotypes.size() > 0 && success)
    {
        SeqItemVector reads;
        SeqItemVector readMates;
        SeqItemVector rcReads;
        SeqItemVector rcReadMates;
        
        // Extract reads from the base BWT
        HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache,
                                          m_parameters.pBaseSSA, m_parameters.kmer, false, &reads, &readMates);

        HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pBaseBWT, m_parameters.pBaseBWTCache,
                                          m_parameters.pBaseSSA, m_parameters.kmer, true, &rcReads, &rcReadMates);

        // Extract reads from the variant BWT
        HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache,
                                          m_parameters.pVariantSSA, m_parameters.kmer, false, &reads, &readMates);

        HapgenUtil::extractHaplotypeReads(inHaplotypes, m_parameters.pVariantBWT, m_parameters.pVariantBWTCache,
                                          m_parameters.pVariantSSA, m_parameters.kmer, true, &rcReads, &rcReadMates);

        //
        // Run Dindel
        //

        double MAP_QUAL = 40.0;
        int BASE_QUAL = 20;
        std::vector<DindelRead> dReads;

        for(size_t i = 0; i < reads.size(); ++i) 
            dReads.push_back(DindelRead(reads[i],std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, true));

        for(size_t i = 0; i < rcReads.size(); ++i)
        {
            rcReads[i].seq.reverseComplement();
            dReads.push_back(DindelRead(rcReads[i],std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, false));
        }
    
        // We need at least 2 haplotypes
        assert(flankingHaplotypes.size() >= 2);
        assert(flankingHaplotypes[0].size() > 0);

        std::string dindelRef = flankingHaplotypes[0]; // First flanking haplotype is of the reference
        
        // The first base of the haplotype sequences on the reference
        int dindelRefStart = 0; 
        // Debug: Separate the non-reference haplotypes from the reference
        try
        {
            StringVector nonReference(flankingHaplotypes.begin()+1, flankingHaplotypes.end());
            std::cout << "Running dindel on " << nonReference.size() << " haplotypes and " << dReads.size() << " reads\n";
            DindelWindow dWindow(nonReference, dindelRef, dindelRefStart, "unknown" );

            DindelRealignParameters dRealignParameters("addSNPMaxSNPs:0");

            DindelRealignWindow dRealignWindow(&dWindow, dReads, dRealignParameters);

            dRealignWindow.run("hmm", m_variantVCFFile);
        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            exit(EXIT_FAILURE);
        }
    }
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
        case BRC_NO_SOLUTION:
            m_stats.numNoSolution += 1;
            break;
    }
    
    return result;
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

//
// GraphCompareAggregateResult
//
GraphCompareAggregateResults::GraphCompareAggregateResults(const std::string& filename) : m_numVariants(0)
{
    //
    m_pWriter = createWriter(filename);

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
}

//
void GraphCompareAggregateResults::printStats() const
{
    m_stats.print();
}
