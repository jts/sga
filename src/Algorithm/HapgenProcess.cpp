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

//
//
//
HapgenProcess::HapgenProcess(const HapgenParameters& params) : m_parameters(params)
{
}

//
HapgenProcess::~HapgenProcess()
{
}

void HapgenProcess::processSite(const std::string& refName, size_t start, size_t end, const std::string& comment)
{
    std::cout << "\n\nProcessing " << refName << " [" << start << " " << end << "] " << comment << "\n";
    std::string startSeq = findAnchorKmer(refName, start, true);
    std::string endSeq = findAnchorKmer(refName, end, false);

    if(startSeq.empty() || endSeq.empty())
        return;

    HaplotypeBuilder builder;
    builder.setTerminals(startSeq, 1, endSeq, 1);
    builder.setIndex(m_parameters.pBWT, m_parameters.pRevBWT);
    builder.setKmerParameters(m_parameters.kmer, m_parameters.kmerThreshold);
    builder.run();
    
    HaplotypeBuilderResult result;
    builder.parseWalks(result);
    std::cout << "Built " << result.haplotypes.size() << " candidate haplotypes\n";

    if(result.haplotypes.size() >= 2)
    {
        std::cout << "Alignment of first two:\n";
        StdAlnTools::globalAlignment(result.haplotypes[0], result.haplotypes[1], true);
    }
}

// Returns the closest kmer to the provided position with occurrence count greater than the passed in threshold
std::string HapgenProcess::findAnchorKmer(const std::string& refName, int64_t position, bool upstream)
{
    const SeqItem& refItem = m_parameters.pRefTable->getRead(refName);
    std::string seq;

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
            return testSeq;
    }

    return "";
}
