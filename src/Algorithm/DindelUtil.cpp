///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
//
// DindelUtil - Wrappers and utility functions
// for the dindel haplotype scoring functions
//
#include "DindelUtil.h"
#include "HapgenUtil.h"

// Run dindel on a pair of samples
DindelReturnCode DindelUtil::runDindelPair(const std::string& normalString, 
                                           const std::string& variantString, 
                                           const GraphCompareParameters& parameters,
                                           std::ostream& baseOut,
                                           std::ostream& variantOut)
{
    StringVector inHaplotypes;
    inHaplotypes.push_back(normalString);
    inHaplotypes.push_back(variantString);

    //
    // First, extract the reads from the normal and variant data sets that match each haplotype
    //

    // Normal reads
    SeqItemVector normalReads;
    SeqItemVector normalReadMates;
    SeqItemVector normalRCReads;
    SeqItemVector normalRCReadMates;
    
    // Set the value to use for extracting reads that potentially match the haplotype
    // Do not use a kmer greater than 41
    size_t KMER_CEILING = 51;
    size_t extractionKmer = parameters.kmer < KMER_CEILING ? parameters.kmer : KMER_CEILING;

    // Reads on the same strand as the haplotype    
    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                      parameters.pBaseSSA, extractionKmer, false, &normalReads, &normalReadMates);

    // Reads on the reverse strand
    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                      parameters.pBaseSSA, extractionKmer, true, &normalRCReads, &normalRCReadMates);

    // Variant reads
    SeqItemVector variantReads;
    SeqItemVector variantReadMates;
    SeqItemVector variantRCReads;
    SeqItemVector variantRCReadMates;

    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                      parameters.pVariantSSA, extractionKmer, false, &variantReads, &variantReadMates);

    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                      parameters.pVariantSSA, extractionKmer, true, &variantRCReads, &variantRCReadMates);


    //
    // Run dindel against the haplotypes for the best alignment position
    // We do this for the normal and variant reads separately
    //
    HapgenAlignment bestAlignment;
    DindelReturnCode retCode = computeBestAlignment(inHaplotypes, variantReadMates, variantRCReadMates, parameters, bestAlignment);

    if(retCode != DRC_OK)
        return retCode;

    // Generate the input haplotypes for dindel
    int FLANKING_SIZE = 0;
    StringVector flankingHaplotypes;
    HapgenUtil::makeFlankingHaplotypes(bestAlignment, parameters.pRefTable, 
                                       FLANKING_SIZE, inHaplotypes, flankingHaplotypes);

    // We need at least 2 haplotypes
    if(flankingHaplotypes.size() < 2)
        return DRC_NO_ALIGNMENT;
    
    // 
    if(flankingHaplotypes[0].size() == 0)
        return DRC_NO_ALIGNMENT;

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

        int dindelRefStart = bestAlignment.position + 1; // VCF coordinates are 1-based
        std::stringstream refName;
        refName << parameters.pRefTable->getRead(bestAlignment.referenceID).id;

        std::string dindelRef = flankingHaplotypes[0]; // First flanking haplotype is of the reference
        StringVector nonReference(flankingHaplotypes.begin()+1, flankingHaplotypes.end());

        try
        {
            DindelWindow dWindow(nonReference, dindelRef, dindelRefStart, refName.str() );

            DindelRealignParameters dRealignParameters("addSNPMaxSNPs:0");
            DindelRealignWindow dRealignWindow(&dWindow, dReads, dRealignParameters);

            std::ostream& currOut = i == 0 ? baseOut : variantOut;
            dRealignWindow.run("hmm", currOut);
        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            return DRC_EXCEPTION;
            //exit(EXIT_FAILURE);
        }
    }

    return DRC_OK;
}

// Compute the best alignment of the haplotype collection to the reference
DindelReturnCode DindelUtil::computeBestAlignment(const StringVector& inHaplotypes, 
                                                  const SeqItemVector& variantMates,
                                                  const SeqItemVector& variantRCMates,
                                                  const GraphCompareParameters& parameters,
                                                  HapgenAlignment& bestAlignment)
{
    size_t MAX_DEPTH = 2000;
    if(variantMates.size() + variantRCMates.size() > MAX_DEPTH)
        return DRC_OVER_DEPTH;

    //
    // Align the haplotypes to the reference genome to generate candidate alignments
    //
    HapgenAlignmentVector candidateAlignments;
    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        HapgenUtil::alignHaplotypeToReference(inHaplotypes[i], parameters.pReferenceBWT, 
                                              parameters.pReferenceSSA, candidateAlignments);
    }

    // Remove duplicate or bad alignment pairs
    HapgenUtil::coalesceAlignments(candidateAlignments);

    if(candidateAlignments.empty())
        return DRC_NO_ALIGNMENT;

    //
    // Score each candidate alignment against the mates of all the variant reads
    //
    int bestCandidate = -1;
    double bestAverageScoreFrac = 0.0f;
    double secondBest = 0.0f;
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        // Compute the average score of the reads' mates to the flanking sequence
        StringVector referenceFlanking;
        HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i], parameters.pRefTable, 
                                           1000, inHaplotypes, referenceFlanking);


        // If valid flanking haplotypes could not be made, skip this alignment
        if(referenceFlanking.empty())
            continue;

        // Realign the mates
        LocalAlignmentResultVector localAlignments = HapgenUtil::alignReadsLocally(referenceFlanking[0], variantMates);
        LocalAlignmentResultVector localAlignmentsRC = HapgenUtil::alignReadsLocally(referenceFlanking[0], variantRCMates);

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

        //printf("Alignment %zu mate-score: %lf\n", i, score);
    }

    //printf("total alignments: %zu best score: %lf\n", candidateAlignments.size(), bestAverageScoreFrac);

    if(bestCandidate == -1)
        return DRC_NO_ALIGNMENT;

    if(bestAverageScoreFrac < 0.9f)
        return DRC_POOR_ALIGNMENT;

    if(bestAverageScoreFrac - secondBest < 0.05f)
        return DRC_AMBIGUOUS_ALIGNMENT;

    bestAlignment = candidateAlignments[bestCandidate];
    return DRC_OK;
}

// Initialize an array of code counts
void DindelUtil::initializeCodeCounts(int counts[DRC_NUM_CODES])
{
    for(int i = 0; i < DRC_NUM_CODES; ++i)
        counts[i] = 0;
}

// Print a report of the dindel return codes
void DindelUtil::printReturnReport(const int counts[DRC_NUM_CODES])
{
    // Sum the total number of haplotype sets processed
    int sum = 0;
    for(int i = 0; i < DRC_NUM_CODES; ++i)
        sum += counts[i];

    printf("Total variants processed: %d\n", sum);
    printf("    number failed due to depth check: %d\n", counts[DRC_OVER_DEPTH]);
    printf("    number failed due to no alignment: %d\n", counts[DRC_NO_ALIGNMENT]);
    printf("    number failed due to poor alignment: %d\n", counts[DRC_POOR_ALIGNMENT]);
    printf("    number failed due to ambiguous alignment: %d\n", counts[DRC_AMBIGUOUS_ALIGNMENT]);
    printf("    number failed due to dindel exception: %d\n", counts[DRC_EXCEPTION]);
    printf("    number passed to dindel: %d\n", counts[DRC_OK]);
}

