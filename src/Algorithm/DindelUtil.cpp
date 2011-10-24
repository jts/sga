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
void DindelUtil::runDindelPair(const std::string& normalString, 
                               const std::string& variantString, 
                               const GraphCompareParameters& parameters,
                               VCFFile& baseVCFFile,
                               VCFFile& variantVCFFile)
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
    
    // Reads on the same strand as the haplotype    
    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                      parameters.pBaseSSA, parameters.kmer, false, &normalReads, &normalReadMates);

    // Reads on the reverse strand
    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                      parameters.pBaseSSA, parameters.kmer, true, &normalRCReads, &normalRCReadMates);

    // Variant reads
    SeqItemVector variantReads;
    SeqItemVector variantReadMates;
    SeqItemVector variantRCReads;
    SeqItemVector variantRCReadMates;

    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                      parameters.pVariantSSA, parameters.kmer, false, &variantReads, &variantReadMates);

    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                      parameters.pVariantSSA, parameters.kmer, true, &variantRCReads, &variantRCReadMates);

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
        std::cout << "Processing " << candidateAlignments[i] << "\n";
        HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i], parameters.pRefTable, 
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
    HapgenUtil::makeFlankingHaplotypes(candidateAlignments[bestCandidate], parameters.pRefTable, 
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
        refName << parameters.pRefTable->getRead(candidateAlignments[bestCandidate].referenceID).id;

        std::string dindelRef = flankingHaplotypes[0]; // First flanking haplotype is of the reference
        StringVector nonReference(flankingHaplotypes.begin()+1, flankingHaplotypes.end());
        std::cout << "Running dindel on " << nonReference.size() << " haplotypes and " << dReads.size() << " reads\n";

        try
        {
            DindelWindow dWindow(nonReference, dindelRef, dindelRefStart, refName.str() );

            DindelRealignParameters dRealignParameters("addSNPMaxSNPs:0");
            DindelRealignWindow dRealignWindow(&dWindow, dReads, dRealignParameters);

            dRealignWindow.run("hmm", i == 0 ? baseVCFFile : variantVCFFile);
        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            exit(EXIT_FAILURE);
        }
    }
}
