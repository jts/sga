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
    size_t KMER_CEILING = 51;
    size_t extractionKmer = parameters.kmer < KMER_CEILING ? parameters.kmer : KMER_CEILING;

    // Reads on the same strand as the haplotype    
    if(!parameters.bReferenceMode)
    {
        HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                          parameters.pBaseSSA, extractionKmer, false, &normalReads, &normalReadMates);

        // Reads on the reverse strand
        HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                          parameters.pBaseSSA, extractionKmer, true, &normalRCReads, &normalRCReadMates);
    }

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

    // If the best alignment is to the reverse strand, flip everything
    if(bestAlignment.isRC)
    {
        // Flip haplotypes
        for(size_t i = 0; i < inHaplotypes.size(); ++i)
            inHaplotypes[i] = reverseComplement(inHaplotypes[i]);

        // Flip alignment strand bit
        // Alignment coordinates do not need to change
        bestAlignment.isRC = false;

        // Swap read vectors
        normalReads.swap(normalRCReads);
        normalReadMates.swap(normalRCReadMates);

        variantReads.swap(variantRCReads);
        variantReadMates.swap(variantRCReadMates);
    }

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
        std::vector< std::vector<DindelReferenceMapping> > dRefMappings(flankingHaplotypes.size());

        // FIXME Set haplotype mapping quality properly
        DindelReferenceMapping dRefMapping(refName.str(), dindelRef, dindelRefStart, 60.0, false);
        for(size_t j = 0; j < flankingHaplotypes.size(); ++j) dRefMappings[j].push_back(dRefMapping);

        try
        {

            DindelWindow dWindow(flankingHaplotypes, dRefMappings);

            DindelRealignParameters dRealignParameters;
            DindelRealignWindow dRealignWindow(&dWindow, dReads, dRealignParameters);

            dRealignWindow.run("hmm", (i==0) ? baseOut : variantOut);

        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            return DRC_EXCEPTION;
        }
    }

    return DRC_OK;
}

// Run dindel on a pair of samples
DindelReturnCode DindelUtil::runDindelPairMatePair(const std::string& normalString,
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
    if(!parameters.bReferenceMode)
    {
        HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                          parameters.pBaseSSA, extractionKmer, false, &normalReads, &normalReadMates);

        // Reads on the reverse strand
        HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                          parameters.pBaseSSA, extractionKmer, true, &normalRCReads, &normalRCReadMates);
    }

    // Variant reads
    SeqItemVector variantReads;
    SeqItemVector variantReadMates;
    SeqItemVector variantRCReads;
    SeqItemVector variantRCReadMates;

    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                      parameters.pVariantSSA, extractionKmer, false, &variantReads, &variantReadMates);

    HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                      parameters.pVariantSSA, extractionKmer, true, &variantRCReads, &variantRCReadMates);

    size_t total_reads = normalReads.size() + normalReadMates.size() + normalRCReads.size() + normalRCReadMates.size();
    total_reads += variantReads.size() + variantReadMates.size() + variantRCReads.size() + variantRCReadMates.size();

    size_t MAX_READS = 1000;

    if(total_reads > MAX_READS)
        return DRC_OVER_DEPTH;

    // Get canidate alignments for the input haplotypes
    HapgenAlignmentVector candidateAlignments;

    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        HapgenAlignmentVector thisCandidateAlignments;
        HapgenUtil::alignHaplotypeToReference(inHaplotypes[i],
                                              parameters.pReferenceBWT,
                                              parameters.pReferenceSSA,
                                              thisCandidateAlignments);

        /*
        std::cout << "Alignments for haplotype " << i << "\n";
        for(size_t j = 0; j < thisCandidateAlignments.size(); ++j)
            std::cout << thisCandidateAlignments[j] << "\n";
        */
        candidateAlignments.insert(candidateAlignments.end(), thisCandidateAlignments.begin(), thisCandidateAlignments.end());
    }
    
    // Remove duplicate or bad alignment pairs
    HapgenUtil::coalesceAlignments(candidateAlignments);
//    std::cout << "Found " << candidateAlignments.size() << " alignments for haplotypes\n";

    // Join each haplotype with flanking sequence from the reference genome for each alignment
    // This function also adds a haplotype (with flanking sequence) for the piece of the reference
    bool success = true;
    int FLANKING_SIZE = 2000;
    StringVector flankingHaplotypes;

    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        // FIXME. Maybe should use only inHaplotypes[i]???
        success = HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i],
                                                     parameters.pRefTable,
                                                     FLANKING_SIZE,
                                                     inHaplotypes,
                                                     flankingHaplotypes);
    
    }


    // Generate the input haplotypes for dindel
    // We need at least 2 haplotypes
    size_t totFlankingHaplotypes = flankingHaplotypes.size();

    if(totFlankingHaplotypes < 2)
        return DRC_NO_ALIGNMENT;

    // Ensure the reference haplotype is a non-empty string
    if(flankingHaplotypes[0].size() == 0)
        return DRC_NO_ALIGNMENT;

    // Make Dindel referenceMappings
    std::vector< std::vector<DindelReferenceMapping> > dindelRefMappings(flankingHaplotypes.size());
    StringVector dindelHaplotypes;
    std::set<DindelReferenceMapping> refMappings;

    //
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        std::string upstream, defined, downstream;
        std::string refName = parameters.pRefTable->getRead(candidateAlignments[i].referenceID).id;

        HapgenUtil::extractReferenceSubstrings(candidateAlignments[i],parameters.pRefTable, FLANKING_SIZE, upstream, defined, downstream);
        std::string refSeq = upstream + defined + downstream;
       
        if(0)
        {
            std::cout << "\n ================================================\n";
            std::cout << "candidateAlignments[" << i << "]" << candidateAlignments[i] << "\n";
            std::cout << "refSeq: " << upstream << " " << defined << " " << downstream << "\n";
            for (size_t k = 0; k < inHaplotypes.size(); ++k)
               std::cout << "        " << std::string(upstream.size(),' ') << " " << inHaplotypes[k] << "\n";

            for (size_t k = 0; k < flankingHaplotypes.size(); ++k)
                std::cout << "        " << " " << flankingHaplotypes[k] << "\n";
        }

        int refStart = candidateAlignments[i].position - int(upstream.size()) + 1;

        // Here the score is used as an estimate of how unique "defined" is in the reference sequence.
        // "defined" is not the reference sequence but a candidate haplotype.
        // It is conservative because the flanking sequence is not used in this estimation.
        DindelReferenceMapping rm(refName, refSeq, refStart, double(candidateAlignments[i].score), candidateAlignments[i].isRC);
        std::set<DindelReferenceMapping>::iterator rmit = refMappings.find(rm);
        if(rmit == refMappings.end())
        {
            refMappings.insert(rm);
        }
        else
        {
            if(rm.referenceAlignmentScore > rmit->referenceAlignmentScore) 
                rmit->referenceAlignmentScore = rm.referenceAlignmentScore;
        }
    }
    
    /*
    std::cout << "REFERENCE MAPPINGS: \n";
    int c = 0;
    for(std::set<DindelReferenceMapping>::const_iterator it = refMappings.begin(); it != refMappings.end(); it++, c++)
        std::cout << c << " " << it->refName << " start: " << it->refStart << " end: " << it->refStart + it->refSeq.size()-1 << " score: " << it->referenceAlignmentScore << "\n";
    */

    for(size_t i = 0; i < flankingHaplotypes.size(); ++i)
    {
        dindelHaplotypes.push_back(flankingHaplotypes[i]);
        dindelRefMappings[i] = std::vector<DindelReferenceMapping>(refMappings.begin(),refMappings.end());
    }

    //
    // Run Dindel
    //
    double MAP_QUAL = 40.0;
    int BASE_QUAL = 20;
    size_t start_i = parameters.bReferenceMode ? 1 : 0;
    for(size_t i = start_i; i <= 1; ++i)
    {
        //std::cout << (i == 0 ? "NORMAL\n" : "VARIANT\n");

        SeqItemVector& fwdReads = (i == 0) ? normalReads : variantReads;
        SeqItemVector& fwdReadMates = (i == 0) ? normalReadMates : variantReadMates;
        SeqItemVector& rcReads = (i == 0) ? normalRCReads : variantRCReads;
        SeqItemVector& rcReadMates = (i == 0) ? normalRCReadMates : variantRCReadMates;

        // Create dindel reads
        // Mates must be at the end of the array.
        std::vector<DindelRead> dReads;
        for(size_t j = 0; j < fwdReads.size(); ++j)
            dReads.push_back(DindelRead(fwdReads[j], std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, true));

        for(size_t j = 0; j < rcReads.size(); ++j)
        {
            rcReads[j].seq.reverseComplement();
            dReads.push_back(DindelRead(rcReads[j], std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, false));
        }

        for(size_t j = 0; j < fwdReadMates.size(); ++j)
        {           
            fwdReadMates[j].seq.reverseComplement();
            dReads.push_back(DindelRead(fwdReadMates[j], std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, true));
        }

        for(size_t j = 0; j < rcReadMates.size(); ++j)
        {       
            //rcReadMates[j].seq.reverseComplement();
            dReads.push_back(DindelRead(rcReadMates[j], std::string("SAMPLE"), MAP_QUAL, BASE_QUAL, false));
        }

        // std::cout << "*******MULTIPLE ALIGNMENT of reads and haplotypes\n";
        // doMultipleReadHaplotypeAlignment(dReads, flankingHaplotypes);

        try
        {
            DindelWindow dWindow(dindelHaplotypes, dindelRefMappings);
            DindelRealignWindow dRealignWindow(&dWindow, dReads, parameters.dindelRealignParameters);
            dRealignWindow.run("hmm", (i==0) ? baseOut : variantOut);
	    //dRealignWindow.run("hmm", std::cout);
        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            exit(DRC_EXCEPTION);
        }

        baseOut.flush();
        variantOut.flush();
    }

    return DRC_OK;
}


void DindelUtil::doMultipleReadHaplotypeAlignment(const std::vector<DindelRead> & dReads,
                                                  const StringVector & haplotypes)
{

    // globally align haplotypes to the first haplotype (arbitrary)
    assert(haplotypes.size()>0);
    for (size_t h = 0; h < haplotypes.size(); ++h)
    {
        std::cout << "ALIGNING EVERYTHING AGAINST HAPLOTYPE " << h << "\n";
        std::vector< MAlignData > maVector;
        const std::string rootSequence = haplotypes[h];
        std::string hid;
        for(size_t j = 0; j < haplotypes.size(); j++)
        {
            MAlignData _ma;
            _ma.position = 0;
            _ma.str = haplotypes[j];

            std::stringstream ss;
            if (j!=h)
                ss << "haplotype-" << j;
            else
                ss << "HAPLOTYPE-" << j;
            _ma.name = ss.str();
            _ma.expandedCigar = StdAlnTools::expandCigar(StdAlnTools::globalAlignmentCigar(haplotypes[j], rootSequence));
            maVector.push_back(_ma);
        }
    

        for(size_t r = 0; r < dReads.size(); ++r)
        {
            MAlignData _ma;
            _ma.position = 0;
            _ma.str = dReads[r].getSequence();

            std::stringstream ss;
            if (r<dReads.size()/2) 
                ss << "read-" << r; 
            else 
                ss << "MATE read-" << r;

            _ma.name = ss.str();
            _ma.expandedCigar = StdAlnTools::expandCigar(StdAlnTools::globalAlignmentCigar(dReads[r].getSequence(), rootSequence));
            maVector.push_back(_ma);
        }

        MultiAlignment MA(rootSequence, maVector, hid);
        std::string consensus = MA.generateConsensus();
        MA.print(100000, &consensus);
    }
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

    if(bestCandidate == -1)
        return DRC_NO_ALIGNMENT;

    /*
    if(bestAverageScoreFrac < 0.9f)
        return DRC_POOR_ALIGNMENT;

    if(bestAverageScoreFrac - secondBest < 0.05f)
        return DRC_AMBIGUOUS_ALIGNMENT;
    */
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

