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
#include "VCFUtil.h"
#include "Profiler.h"

// Run dindel on a pair of samples
DindelReturnCode DindelUtil::runDindelPairMatePair(const std::string& id,
                                                   const StringVector& base_haplotypes,
                                                   const StringVector& variant_haplotypes,
                                                   const GraphCompareParameters& parameters,
                                                   std::ostream& baseOut,
                                                   std::ostream& variantOut)
{
    PROFILE_FUNC("runDindelPairMatePair")

    StringVector inHaplotypes;
    inHaplotypes.insert(inHaplotypes.end(), base_haplotypes.begin(), base_haplotypes.end());
    inHaplotypes.insert(inHaplotypes.end(), variant_haplotypes.begin(), variant_haplotypes.end());

    //
    // First, extract the reads from the normal and variant data sets that match each haplotype
    //
    size_t MAX_READS = 40000000000;

    // Get canidate alignments for the input haplotypes
    HapgenAlignmentVector candidateAlignments;

    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        HapgenAlignmentVector thisCandidateAlignments;
        HapgenUtil::alignHaplotypeToReferenceKmer(inHaplotypes[i],
                                                  parameters.pReferenceBWT,
                                                  parameters.pReferenceSSA,
                                                  parameters.pRefTable,
                                                  thisCandidateAlignments);

        candidateAlignments.insert(candidateAlignments.end(), thisCandidateAlignments.begin(), thisCandidateAlignments.end());
    }
    
    // Remove duplicate or bad alignment pairs
    HapgenUtil::coalesceAlignments(candidateAlignments);
    
    size_t MAX_ALIGNMENTS = 10;
    printf("Found %zu alignments\n", candidateAlignments.size());
    if(candidateAlignments.size() > MAX_ALIGNMENTS)
        return DRC_AMBIGUOUS_ALIGNMENT;

    // Join each haplotype with flanking sequence from the reference genome for each alignment
    // This function also adds a haplotype (with flanking sequence) for the piece of the reference
    bool success = true;
    int FLANKING_SIZE = 1000;
    StringVector flankingHaplotypes;

    // This vector contains the internal portion of the haplotypes, without the flanking sequence
    // It is used to extract reads
    StringVector candidateHaplotypes;
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        // FIXME. Maybe should use only inHaplotypes[i]???
        success = HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i],
                                                     parameters.pRefTable,
                                                     FLANKING_SIZE,
                                                     inHaplotypes,
                                                     flankingHaplotypes,
                                                     candidateHaplotypes);
    
    }

    // Normal reads
    SeqItemVector normalReads;
    SeqItemVector normalReadMates;
    SeqItemVector normalRCReads;
    SeqItemVector normalRCReadMates;

    // Remove non-unique candidate haplotypes
    std::sort(candidateHaplotypes.begin(), candidateHaplotypes.end());
    StringVector::iterator haplotype_iterator = std::unique(candidateHaplotypes.begin(), candidateHaplotypes.end());
    candidateHaplotypes.resize(haplotype_iterator - candidateHaplotypes.begin());

    // Set the value to use for extracting reads that potentially match the haplotype
    // Do not use a kmer for extraction greater than this value
    size_t KMER_CEILING = 41;
    size_t extractionKmer = parameters.kmer < KMER_CEILING ? parameters.kmer : KMER_CEILING;
    
    bool extractOK = true;
    if(!parameters.bReferenceMode)
    {
        // Reads on the same strand as the haplotype
        extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                                              parameters.pBaseSSA, extractionKmer, false, MAX_READS, &normalReads, &normalReadMates);

        if(!extractOK)
            return DRC_OVER_DEPTH;

        // Reads on the reverse strand
        extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                                              parameters.pBaseSSA, extractionKmer, true, MAX_READS, &normalRCReads, &normalRCReadMates);

        if(!extractOK)
            return DRC_OVER_DEPTH;
    }

    // Variant reads
    SeqItemVector variantReads;
    SeqItemVector variantReadMates;
    SeqItemVector variantRCReads;
    SeqItemVector variantRCReadMates;

    extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                                  parameters.pVariantSSA, extractionKmer, false, MAX_READS, &variantReads, &variantReadMates);

    if(!extractOK)
        return DRC_OVER_DEPTH;

    extractOK = HapgenUtil::extractHaplotypeReads(candidateHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                                  parameters.pVariantSSA, extractionKmer, true, MAX_READS, &variantRCReads, &variantRCReadMates);

    if(!extractOK)
        return DRC_OVER_DEPTH;

    size_t total_reads = normalReads.size() + normalReadMates.size() + normalRCReads.size() + normalRCReadMates.size();
    total_reads += variantReads.size() + variantReadMates.size() + variantRCReads.size() + variantRCReadMates.size();

    if(total_reads > MAX_READS)
        return DRC_OVER_DEPTH;    

    printf("Passing to dindel %zu haplotypes, %zu reads\n", candidateAlignments.size(), total_reads);

    // Generate the input haplotypes for dindel
    // We need at least 2 haplotypes
    size_t totFlankingHaplotypes = flankingHaplotypes.size();

    if(totFlankingHaplotypes < 2)
        return DRC_NO_ALIGNMENT;

    // Ensure the reference haplotype is a non-empty string
    if(flankingHaplotypes[0].size() == 0)
        return DRC_NO_ALIGNMENT;


    // Make Dindel referenceMappings
    // std::vector< std::vector<DindelReferenceMapping> > dindelRefMappings(flankingHaplotypes.size());
    StringVector dindelHaplotypes;
    std::set<DindelReferenceMapping> refMappings;

    //
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        std::string upstream, defined, downstream;
        std::string refName = parameters.pRefTable->getRead(candidateAlignments[i].referenceID).id;

        HapgenUtil::extractReferenceSubstrings(candidateAlignments[i],parameters.pRefTable, FLANKING_SIZE, upstream, defined, downstream);
        std::string refSeq = upstream + defined + downstream;
     
        int refStart = candidateAlignments[i].position - int(upstream.size()) + 1;

        // Here the score is used as an estimate of how unique "defined" is in the reference sequence.
        // "defined" is not the reference sequence but a candidate haplotype.
        // It is conservative because the flanking sequence is not used in this estimation.
        // DindelReferenceMapping rm(refName, refSeq, refStart, double(candidateAlignments[i].score), candidateAlignments[i].isRC);
    	DindelReferenceMapping rm(refName, refSeq, refStart, double(candidateAlignments[i].score+2*FLANKING_SIZE), candidateAlignments[i].isRC);
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
    
     
    std::cout << "REFERENCE MAPPINGS: \n";
    int c = 0;
    for(std::set<DindelReferenceMapping>::const_iterator it = refMappings.begin(); it != refMappings.end(); it++, c++)
        std::cout << c << " " << it->refName << " start: " << it->refStart << " end: " << it->refStart + it->refSeq.size()-1 << " score: " << it->referenceAlignmentScore << "\n";
    
    for(size_t i = 0; i < flankingHaplotypes.size(); ++i)
    {
        dindelHaplotypes.push_back(flankingHaplotypes[i]);
        //dindelRefMappings[i] = std::vector<DindelReferenceMapping>(refMappings.begin(),refMappings.end());
    }
    std::vector<DindelReferenceMapping> dRefMappings(refMappings.begin(),refMappings.end());
    DindelWindow dWindow(dindelHaplotypes, dRefMappings);

    if (0)
    {
        for (size_t i = 0; i < dindelHaplotypes.size(); i++ )
        {
            std::cout << ">HAPLOTYPE_" << i << "\n";
            std::cout << dindelHaplotypes[i] << "\n";
        }
    }

    //
    // Run Dindel
    //
    double MAP_QUAL = 40.0;
    int BASE_QUAL = 20;
    size_t start_i = parameters.bReferenceMode ? 1 : 0;

    DindelRealignWindowResult *pThisResult, *pPreviousResult = NULL;

    for(size_t i = start_i; i <= 1; ++i)
    {

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

        if (parameters.dindelRealignParameters.realignMatePairs)
        {
            std::cout << "Adding read mates.\n";
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
        }

//        std::cout << "*******MULTIPLE ALIGNMENT of reads and haplotypes\n";
//        doMultipleReadHaplotypeAlignment(dReads, flankingHaplotypes);

        pThisResult = new DindelRealignWindowResult();

        try
        {
            DindelRealignWindow dRealignWindow(&dWindow, dReads, parameters.dindelRealignParameters);
            dRealignWindow.run("hmm", (i==0) ? baseOut : variantOut, id, pThisResult, pPreviousResult);
	    //dRealignWindow.run("hmm", std::cout);
        }
        catch(std::string e)
        {
            std::cerr << "Dindel Exception: " << e << "\n";
            exit(DRC_EXCEPTION);
        }

        baseOut.flush();
        variantOut.flush();

        if(i == 0)
            pPreviousResult = pThisResult;
    }

    delete pThisResult;
    delete pPreviousResult;

    return DRC_OK;
}

// Run dindel on a pair of samples
DindelReturnCode DindelUtil::runNaiveCaller(const std::string& normalString,
                                            const std::string& variantString,
                                            const GraphCompareParameters& parameters,
                                            std::ostream& baseOut,
                                            std::ostream& variantOut)
{
    int verbose = 0;
    StringVector inHaplotypes;
    inHaplotypes.push_back(normalString);
    inHaplotypes.push_back(variantString);

    //
    // Find alignment of haplotypes to reference
    //
    // Find the best alignment location of the haplotypes
    HapgenAlignmentVector candidateAlignments;
    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        HapgenAlignmentVector thisCandidateAlignments;
        HapgenUtil::alignHaplotypeToReferenceBWASW(inHaplotypes[i],
                                                   parameters.pReferenceBWT,
                                                   parameters.pReferenceSSA,
                                                   candidateAlignments);
    }
    
    // Remove duplicate or bad alignment pairs
    HapgenUtil::coalesceAlignments(candidateAlignments);

    if(candidateAlignments.empty())
        return DRC_NO_ALIGNMENT;

    HapgenAlignment bestAlignment;
    int bestScore = -1;
    for(size_t i = 0; i < candidateAlignments.size(); ++i)
    {
        if(candidateAlignments[i].score > bestScore)
        {
            bestScore = candidateAlignments[i].score;
            bestAlignment = candidateAlignments[i];
        }
    }

    // If the best alignment is wrt reverse strand, flip the haplotypes
    if(bestAlignment.isRC)
    {
        for(size_t i = 0; i < inHaplotypes.size(); ++i)
            inHaplotypes[i] = reverseComplement(inHaplotypes[i]);
    }

    //
    // Extract reads matching haplotypes
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
    size_t MAX_READS = 1000;

    bool extractOK = true;
    // Reads on the same strand as the haplotype
    if(!parameters.bReferenceMode)
    {
        extractOK = HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                                      parameters.pBaseSSA, extractionKmer, false, MAX_READS, &normalReads, &normalReadMates);

        if(!extractOK)
            return DRC_OVER_DEPTH;

        // Reads on the reverse strand
        extractOK = HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pBaseBWT, parameters.pBaseBWTCache,
                                                      parameters.pBaseSSA, extractionKmer, true, MAX_READS, &normalRCReads, &normalRCReadMates);

        if(!extractOK)
            return DRC_OVER_DEPTH;
    }

    // Variant reads
    SeqItemVector variantReads;
    SeqItemVector variantReadMates;
    SeqItemVector variantRCReads;
    SeqItemVector variantRCReadMates;

    extractOK = HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                                  parameters.pVariantSSA, extractionKmer, false, MAX_READS, &variantReads, &variantReadMates);

    if(!extractOK)
        return DRC_OVER_DEPTH;


    extractOK = HapgenUtil::extractHaplotypeReads(inHaplotypes, parameters.pVariantBWT, parameters.pVariantBWTCache,
                                                  parameters.pVariantSSA, extractionKmer, true, MAX_READS, &variantRCReads, &variantRCReadMates);

    if(!extractOK)
        return DRC_OVER_DEPTH;

    //
    // Construct a multiple alignment between the reference, the assembled haplotypes and all the reads
    // 
    int FLANKING = 0;
    StringVector finalHaplotypes;

    // Add the reference haplotype
    finalHaplotypes.push_back(HapgenUtil::extractReference(bestAlignment, parameters.pRefTable, FLANKING));
    assert(FLANKING == 0);
    finalHaplotypes.insert(finalHaplotypes.end(), inHaplotypes.begin(), inHaplotypes.end());

    // Construct a vector of all sequences
    StringVector allSequences;

    // Add the haplotypes
    allSequences.insert(allSequences.end(), finalHaplotypes.begin(), finalHaplotypes.end());

    size_t first_normal_read_idx = allSequences.size();

    for(size_t i = 0; i < normalReads.size(); ++i)
        allSequences.push_back(normalReads[i].seq.toString());
    for(size_t i = 0; i < normalRCReads.size(); ++i)
        allSequences.push_back(reverseComplement(normalRCReads[i].seq.toString()));
    
    // Save the position of the first read for the variant sequences
    size_t first_variant_read_idx = allSequences.size();

    for(size_t i = 0; i < variantReads.size(); ++i)
        allSequences.push_back(variantReads[i].seq.toString());
    for(size_t i = 0; i < variantRCReads.size(); ++i)
        allSequences.push_back(reverseComplement(variantRCReads[i].seq.toString()));

    const std::string& reference_string = allSequences[0];
    MAlignDataVector madVector;
    for(size_t i = 1; i < allSequences.size(); ++i)
    {
        // The alignment is (slightly counterintuitively) with respect to the 
        // base string as the query so that all the cigar strings used
        // are based on edit operations to the base.
        const std::string& seq = allSequences[i];
        LocalAlignmentResult result = StdAlnTools::localAlignment(seq, reference_string);
        
        // Initialize the multiple alignment data
        MAlignData maData;
        maData.position = result.queryStartPosition;
        // If the non-base sequence overhangs the base, clip it
        if(result.targetStartPosition > 0)
            maData.str = seq.substr(result.targetStartPosition);
        else
            maData.str = seq;

        // Pad out the cigar with reference skip characters
        maData.expandedCigar = StdAlnTools::expandCigar(result.cigar);

        std::stringstream name;
        if(i < first_normal_read_idx)
            name << "haplotype-" << i;
        else if(i < first_variant_read_idx)
            name << "normal-" << i - first_normal_read_idx;
        else
            name << "variant-" << i - first_variant_read_idx;

        maData.name = name.str();
        madVector.push_back(maData);
    }

    MultiAlignment ma(reference_string, madVector, "reference");
    if(verbose)
        ma.print(10000, NULL, false, true);

    // Parse variants between the reference and the variant haplotype (idx 2)
    size_t refRow = 0;
    size_t baseRow = 1;
    size_t varRow = 2;
    size_t numCols = ma.getNumColumns();

    std::string refName = parameters.pRefTable->getRead(bestAlignment.referenceID).id;
    int refPosition = bestAlignment.position;

    int eventStart = -1;
    int eventEnd = -1;
    bool isIndel = false;
    bool hasIndel = false;
    for(size_t i = 0; i < numCols; ++i)
    {
        char refSymbol = ma.getSymbol(refRow, i);
        char varSymbol = ma.getSymbol(varRow, i);
        bool isVariant = (varSymbol != refSymbol);
        
        // Process variants indicated by this portion of the alignment
        if(isVariant)
        {
            if(eventStart == -1)
            {
                // Not currently in a variant event, start a new one
                eventStart = i;
                eventEnd = i;
            }
            else
            {
                // Extend the event
                assert(eventEnd != -1);
                eventEnd += 1;
            }

            isIndel = isIndel || varSymbol == '-' || refSymbol == '-';
        }
        else
        {
            // Check if this is the end of a variant
            if(eventStart != -1)
            {
                // Extract the substrings of the three sequences describing
                // this event.
                // If the event is an indel, we move the start point
                // of the event backwards until all 3 columns are not padding characters.
                // This is to avoid the following situation:
                // AATATTT--CGT ref
                // AATATTTT-CGT base
                // AATATTTTTCGT var
                // or 
                // TTAGGTTTTTTT ref 
                // TTAGG-TTTTTT base 
                // TTAGG--TTTTT var
                //
                // In the first case, it is a 1 base insertion wrt base but a 2 base insertion wrt to the ref.
                // So we need to move the event to the first column with all Ts to properly capture what is happening.
                // In the second case it is a 1 base deletion wrt base and a 2 base deletion wrt to the ref. 
                // Again, we need to move back to the last full column.
                if(isIndel)
                {
                    // It is impossible to extract a reference string if the indel is at the beginning
                    // or end of th multiple alignment so we just fail out here
                    if(eventStart == 0 || eventEnd == (int)numCols - 1)
                        return DRC_EXCEPTION;
                    
                    // Move the start of the event backwards until not in a gap
                    while(eventStart >= 0)
                    {
                        char rc = ma.getSymbol(refRow, eventStart);
                        char bc = ma.getSymbol(baseRow, eventStart);
                        char vc = ma.getSymbol(varRow, eventStart);
                        if(rc == '-' || bc == '-' || vc == '-')
                            eventStart -= 1;
                        else
                            break;
                    }
                    assert(eventStart >= 0);
                }

                size_t eventLength = eventEnd - eventStart + 1;
                std::string refString = StdAlnTools::unpad(ma.getPaddedSubstr(refRow, eventStart, eventLength));
                std::string varString = StdAlnTools::unpad(ma.getPaddedSubstr(varRow, eventStart, eventLength));

                hasIndel = hasIndel || isIndel;
                assert(!refString.empty());
                //assert(!baseString.empty());
                assert(!varString.empty());

                // Calculate how many of the normal/variant reads match the variant
                int normal_matches_variant = 0;
                int normal_matches_ref = 0;
                int variant_matches_variant = 0;
                int variant_matches_ref = 0;

                for(size_t idx = first_normal_read_idx; idx < allSequences.size(); ++idx)
                {
                    bool is_normal_read = idx < first_variant_read_idx;

                    // Get the read sequence at the event
                    std::string readString = StdAlnTools::unpad(ma.getPaddedSubstr(idx, eventStart, eventLength));
                    if(readString == varString)
                    {
                        if(is_normal_read)
                            normal_matches_variant += 1;
                        else
                            variant_matches_variant += 1;
                    }

                    if(readString == refString)
                    {
                        if(is_normal_read)
                            normal_matches_ref += 1;
                        else
                            variant_matches_ref += 1;
                    }

                }
                
                if(verbose)
                {
                    printf("Ref start: %d col: %d eventStart: %d eventEnd: %d\n", (int)refPosition, (int)i, eventStart, eventEnd);
                    std::cout << "RefString " << refString << "\n";
                    std::cout << "varString " << varString << "\n";
                    std::cout << "Normal match: " << normal_matches_variant << "\n";
                    std::cout << "Variant match: " << variant_matches_variant << "\n";
                }

                // Make the VCF record
                int refBaseOffset = ma.getBaseIdx(refRow, eventStart);
                
                // Generate VCF Record
                VCFRecord varRecord;
                varRecord.refName = refName;
                varRecord.refPosition = refPosition + refBaseOffset + 1; // VCF uses 1-based coordinates

                varRecord.refStr = refString;
                varRecord.varStr = varString;
                varRecord.setQuality(10.0f * variant_matches_variant); // awful
                varRecord.setPassStr(variant_matches_variant >= 1 ? "PASS" : "NoCall");
                varRecord.addComment("DP", variant_matches_variant);
                varRecord.addComment("DR", variant_matches_ref);
                double af = (double)variant_matches_variant / (variant_matches_variant + variant_matches_ref);
                varRecord.addComment("AF", af);

                // Generate VCF Record
                VCFRecord baseRecord;
                baseRecord.refName = refName;
                baseRecord.refPosition = refPosition + refBaseOffset + 1; // VCF uses 1-based coordinates

                baseRecord.refStr = refString;
                baseRecord.varStr = varString;
                baseRecord.setQuality(10.0f * normal_matches_variant);
                baseRecord.setPassStr(normal_matches_variant >= 1 ? "PASS" : "NoCall");
                baseRecord.addComment("DP", normal_matches_variant);
                baseRecord.addComment("DR", normal_matches_ref);

                std::cout << baseRecord << "\n";
                std::cout << varRecord << "\n";

                variantOut << varRecord << "\n";
                baseOut << baseRecord << "\n";

                // Reset state
                eventStart = -1;
                eventEnd = -1;
                isIndel = false;                
            }
        }
    }

    return DRC_OK;
}

//
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
                ss << "read-" << r << "("  << dReads[r].getID() << ")"; 
            else 
                ss << "MATE read-" << r;

            _ma.name = ss.str();
            _ma.expandedCigar = StdAlnTools::expandCigar(StdAlnTools::globalAlignmentCigar(dReads[r].getSequence(), rootSequence));
            maVector.push_back(_ma);
        }

        MultiAlignment MA(rootSequence, maVector, hid);
        //std::string consensus = MA.generateConsensus();
        MA.print(100000, NULL, true, true);
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
        HapgenUtil::alignHaplotypeToReferenceBWASW(inHaplotypes[i], parameters.pReferenceBWT, 
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
        StringVector referenceHaplotypes;
        HapgenUtil::makeFlankingHaplotypes(candidateAlignments[i], parameters.pRefTable, 
                                           1000, inHaplotypes, referenceFlanking, referenceHaplotypes);

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

