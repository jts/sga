///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HapgenProcess - Generate candidate haplotypes
// from an assembly graph for a stream of sites
//
#include <algorithm>
#include "HapgenProcess.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"
#include "HaplotypeBuilder.h"
#include "MultiAlignment.h"

//
//
//
HapgenProcess::HapgenProcess(const HapgenParameters& params) : m_parameters(params), m_vcfFile(params.vcfOutfile,"w")
{

}

//
HapgenProcess::~HapgenProcess()
{

}

void HapgenProcess::processSite(const std::string& refName, size_t start, size_t end, const std::string& comment)
{
    assert(false);
    (void)refName;
    (void)start;
    (void)end;
    (void)comment;
#if 0

    if(m_parameters.verbose > 0)
        std::cout << "\nProcessing " << refName << " [" << start << " " << end << "] " << comment << "\n";

    AnchorSequence startAnchor = findAnchorKmer(refName, start, true);
    AnchorSequence endAnchor = findAnchorKmer(refName, end, false);

    if(startAnchor.sequence.empty() || endAnchor.sequence.empty())
    {
        if(m_parameters.verbose > 0)
            std::cout << "Could not anchor to reference\n";
        return;
    }

    if(m_parameters.verbose > 0)
    {
        std::cout << "Left anchor depth: " << startAnchor.count << "\n";
        std::cout << "Right anchor depth: " << endAnchor.count << "\n";
    }

    HaplotypeBuilder builder;
    builder.setTerminals(startAnchor, endAnchor);
    builder.setIndex(m_parameters.pBWT, m_parameters.pRevBWT);
    builder.setKmerParameters(m_parameters.kmer, m_parameters.kmerThreshold);
    builder.run();
    
    HaplotypeBuilderResult result;
    builder.parseWalks(result);
    
    if(m_parameters.verbose > 0)
        std::cout << "Built " << result.haplotypes.size() << " candidate haplotypes\n";

    // Extract the reference sequence spanned by the anchors
    const SeqItem& refItem = m_parameters.pRefTable->getRead(refName);
    size_t refStart = startAnchor.position;
    size_t refEnd = endAnchor.position + m_parameters.kmer;
    std::string refSubstring = refItem.seq.substr(refStart, refEnd - refStart);

    if(result.haplotypes.size() >= 2 && m_parameters.verbose > 0)
    {
        SeqItemVector seqVector;
        std::stringstream rssName;
        rssName << refName << ":" << refStart << "-" << refEnd;
        SeqItem rsi = { rssName.str(), refSubstring };
        seqVector.push_back(rsi);

        for(size_t i = 0; i < result.haplotypes.size(); ++i)
        {
            std::stringstream namer;
            namer << "haplotype-" << i;
            SeqItem hsi = { namer.str(), result.haplotypes[i] };
            seqVector.push_back(hsi);
        }

        MultiAlignment haplotypeAlignment = MultiAlignmentTools::alignSequencesGlobal(seqVector);
        haplotypeAlignment.print();
    }


    // Testing code for extracting reads with k-mer matches to the haplotypes.
    SeqItemVector reads;
    SeqItemVector readMates;
    SeqItemVector rcReads;
    SeqItemVector rcReadMates;

    extractHaplotypeReads(result.haplotypes, false, &reads, &readMates);
    extractHaplotypeReads(result.haplotypes, true, &rcReads, &rcReadMates);

    // FIXME
    double mappingQual = 40.0;
    int fixedBaseQual = 20;

    // make a vector of Dindel reads
    std::vector<DindelRead> dReads;
    for(size_t i = 0; i < reads.size(); ++i) dReads.push_back(DindelRead(reads[i],std::string("SAMPLE"), mappingQual, fixedBaseQual, true ));
    for(size_t i = 0; i < rcReads.size(); ++i)
    {
        rcReads[i].seq.reverseComplement();
        dReads.push_back(DindelRead(rcReads[i],std::string("SAMPLE"), mappingQual, fixedBaseQual, false ));
    }
    // create Dindel window

    std::string dindelRef = refSubstring;
    int dindelRefStart = int(refStart);
    std::vector<DindelReferenceMapping> dRefMappings;

    // FIXME Set haplotype mapping quality properly
    DindelReferenceMapping dRefMapping(refName, dindelRef, dindelRefStart, 100.0, false);
    dRefMappings.push_back(dRefMapping);

    DindelWindow dWindow(result.haplotypes, dRefMappings);

    DindelRealignParameters dRealignParameters;
    DindelRealignWindow dRealignWindow(&dWindow, dReads, dRealignParameters);

    assert(false);
//    dRealignWindow.run("hmm", m_vcfFile.getOutputStream());


    if(m_parameters.verbose > 0)
    {
        // Align mates to neighborhood
        size_t window = 1000;
        size_t ns = refStart > window ? refStart - window : 0;
        std::string neighborhood = refItem.seq.substr(ns, 2*window + refSubstring.length());

        for(size_t i = 0; i < readMates.size(); ++i)
        {
            std::cout << "Aligning mate " << readMates[i].id << "\n";
            LocalAlignmentResult localResult = StdAlnTools::localAlignment(neighborhood, readMates[i].seq.toString());
            LocalAlignmentResult rcLocalResult = StdAlnTools::localAlignment(neighborhood, reverseComplement(readMates[i].seq.toString()));
            std::cout << "  result(ss): " << localResult << "\n";
            std::cout << "  result(rc): " << rcLocalResult << "\n";
        }
    }
#endif
}

// Returns the closest kmer to the provided position with occurrence count greater than the passed in threshold
AnchorSequence HapgenProcess::findAnchorKmer(const std::string& refName, int64_t position, bool upstream)
{
    const SeqItem& refItem = m_parameters.pRefTable->getRead(refName);
    AnchorSequence anchor;

    int64_t stride = upstream ? -1 : 1;
    int MAX_DISTANCE = 100;
    int64_t stop = upstream ? position - MAX_DISTANCE : position + MAX_DISTANCE;

    // Cap the travel distance to avoid out of bounds
    if(stop < 0)
        stop = 0;
    if(stop > (int64_t)(refItem.seq.length() - m_parameters.kmer))
        stop = refItem.seq.length() - m_parameters.kmer;

    for(; position != stop; position += stride)
    {
        std::string testSeq = refItem.seq.substr(position, m_parameters.kmer);
        std::transform(testSeq.begin(), testSeq.end(), testSeq.begin(), ::toupper);
        if(testSeq.find_first_of('N') != std::string::npos)
            continue;

        size_t count = BWTAlgorithms::countSequenceOccurrencesWithCache(testSeq, m_parameters.pBWT, m_parameters.pBWTCache);
        if(count > m_parameters.kmerThreshold)
        {
            anchor.sequence = testSeq;
            anchor.count = count;
            anchor.position = position;
            return anchor;
        }
    }

    anchor.sequence = "";
    anchor.count = -1;
    return anchor;
}
