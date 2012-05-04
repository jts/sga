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

// Extract the reads from the FM-index that share a kmer with any given haplotype
void HapgenProcess::extractHaplotypeReads(const StringVector& haplotypes, bool doReverseComp, 
                                         SeqItemVector* pOutReads, SeqItemVector* pOutMates) const
{
    WARN_ONCE("HapgenProcess::extractHaplotypeReads is depracated")
    assert(false);
    // Extract the set of reads that have at least one kmer shared with these haplotypes
    // This is a bit of a roundabout procedure with a few steps:
    // 1) extract all the kmers in the haplotypes
    // 2) find the intervals for the kmers in the fm-index
    // 3) compute the set of read indices of the reads from the intervals (using the sampled suffix array)
    // 4) finally, extract the read sequences from the index

    // Make a set of kmers from the haplotypes
    std::set<std::string> kmerSet;
    for(size_t i = 0; i < haplotypes.size(); ++i)
    {
        const std::string& h = haplotypes[i];
        for(size_t j = 0; j < h.size() - m_parameters.kmer + 1; ++j)
        {
            std::string ks = h.substr(j, m_parameters.kmer);
            if(doReverseComp)
                ks = reverseComplement(ks);
            kmerSet.insert(ks);
        }
    }

    // Compute the set of reads ids 
    std::set<int64_t> readIndices;
    for(std::set<std::string>::const_iterator iter = kmerSet.begin(); iter != kmerSet.end(); ++iter)
    {
        BWTInterval interval = BWTAlgorithms::findIntervalWithCache(m_parameters.pBWT, m_parameters.pBWTCache, *iter);
        for(int64_t j = interval.lower; j <= interval.upper; ++j)
        {
            // Get index from sampled suffix array
            SAElem elem = m_parameters.pSSA->calcSA(j, m_parameters.pBWT);
            readIndices.insert(elem.getID());
        }
    }

    for(std::set<int64_t>::const_iterator iter = readIndices.begin(); iter != readIndices.end(); ++iter)
    {
        int64_t idx = *iter;
        
        // Extract the read
        std::stringstream namer;
        namer << "idx-" << idx;
        SeqItem item;
        item.id = namer.str();
        item.seq = BWTAlgorithms::extractString(m_parameters.pBWT, idx);
        pOutReads->push_back(item);

        // Optionally extract its mate
        // If the index is constructed properly, 
        // paired reads are in adjacent indices with the
        // first read at even indices
        if(pOutMates != NULL)
        {
            int64_t mateIdx = idx;
            if(idx % 2 == 0)
                mateIdx += 1;
            else
                mateIdx -= 1;

            std::stringstream mateName;
            mateName << "idx-" << mateIdx;
            SeqItem mateItem;
            mateItem.id = mateName.str();
            mateItem.seq = BWTAlgorithms::extractString(m_parameters.pBWT, idx);
            pOutMates->push_back(mateItem);
        }
    }
}

