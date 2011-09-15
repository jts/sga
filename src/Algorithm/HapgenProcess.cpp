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
HapgenProcess::HapgenProcess(const HapgenParameters& params) : m_parameters(params)
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

    std::pair<std::string, int> startAnchor = findAnchorKmer(refName, start, true);
    std::pair<std::string, int> endAnchor = findAnchorKmer(refName, end, false);

    if(startAnchor.first.empty() || endAnchor.first.empty())
    {
        if(m_parameters.verbose > 0)
            std::cout << "Could not anchor to reference\n";
        return;
    }

    if(m_parameters.verbose > 0)
    {
        std::cout << "Left anchor depth: " << startAnchor.second << "\n";
        std::cout << "Right anchor depth: " << endAnchor.second << "\n";
    }

    HaplotypeBuilder builder;
    builder.setTerminals(startAnchor.first, startAnchor.second, endAnchor.first, endAnchor.second);
    builder.setIndex(m_parameters.pBWT, m_parameters.pRevBWT);
    builder.setKmerParameters(m_parameters.kmer, m_parameters.kmerThreshold);
    builder.run();
    
    HaplotypeBuilderResult result;
    builder.parseWalks(result);
    
    if(m_parameters.verbose > 0)
        std::cout << "Built " << result.haplotypes.size() << " candidate haplotypes\n";

    if(result.haplotypes.size() >= 2 && m_parameters.verbose > 0)
    {
        MultiAlignment haplotypeAlignment = MultiAlignmentTools::alignSequences(result.haplotypes);
        haplotypeAlignment.print();
    }


    // Testing code for extracting reads with k-mer matches to the haplotypes.
    /*
    StringVector reads;
    StringVector rcReads;
    extractHaplotypeReads(result.haplotypes, reads, false);
    extractHaplotypeReads(result.haplotypes, rcReads, true);

    if(m_parameters.verbose > 0)
    {
        printf("Found %zu reads matching a kmer with a haplotype\n", reads.size());
        printf("Found %zu reads matching a reverse-complement kmer with a haplotype\n", rcReads.size());
        
        if(!result.haplotypes.empty())
        {
            std::cout << "Printing pairwise alignments of forward reads to haplotype 0\n";
            for(size_t i = 0; i < reads.size(); ++i)
            {
                StdAlnTools::globalAlignment(result.haplotypes[0], reads[i], true);
            }
        }            
    }
    */
}

// Returns the closest kmer to the provided position with occurrence count greater than the passed in threshold
std::pair<std::string, int> HapgenProcess::findAnchorKmer(const std::string& refName, int64_t position, bool upstream)
{
    const SeqItem& refItem = m_parameters.pRefTable->getRead(refName);
    std::pair<std::string, int> anchor;

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
            anchor.first = testSeq;
            anchor.second = count;
            return anchor;
        }
    }

    anchor.first = "";
    anchor.second = -1;
    return anchor;
}

// Extract the reads from the FM-index that share a kmer with any given haplotype
void HapgenProcess::extractHaplotypeReads(const StringVector& haplotypes, StringVector& reads, bool doReverseComp) const
{
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
        reads.push_back(BWTAlgorithms::extractString(m_parameters.pBWT, *iter));
}

