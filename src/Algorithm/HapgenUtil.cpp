///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// HapgenUtil - Utility functions used for 
// working with the haplotype generation
// modules and output
//

#include "HapgenUtil.h"
#include "LRAlignment.h"
#include "Interval.h"
#include "Profiler.h"

// Align the haplotype to the reference genome represented by the BWT/SSA pair
void HapgenUtil::alignHaplotypeToReference(const std::string& haplotype,
                                           const BWT* pReferenceBWT,
                                           const SampledSuffixArray* pReferenceSSA,
                                           HapgenAlignmentVector& outAlignments)
{
    PROFILE_FUNC("HapgenUtil::alignHaplotypesToReference")
    LRAlignment::LRParams params;

    params.zBest = 5;

    for(size_t i = 0; i <= 1; ++i)
    {
        LRAlignment::LRHitVector hits;
        std::string query = (i == 0) ? haplotype : reverseComplement(haplotype);
        LRAlignment::bwaswAlignment(query, pReferenceBWT, pReferenceSSA, params, hits);

        // Convert the hits into alignments
        for(size_t j = 0; j < hits.size(); ++j)
        {
            int q_alignment_length = hits[j].q_end - hits[j].q_start;

            // Skip non-complete alignments
            if((int)haplotype.length() == q_alignment_length)
            {
                HapgenAlignment aln(hits[j].targetID, hits[j].t_start, hits[j].length, hits[j].G, i == 1);
                outAlignments.push_back(aln);
            }
        }
    }
}

// Coalesce a set of alignments into distinct locations
void HapgenUtil::coalesceAlignments(HapgenAlignmentVector& alignments)
{
    if(alignments.empty())
        return;

    // Sort the alignments by reference id, then position
    std::sort(alignments.begin(), alignments.end());

    HapgenAlignmentVector outAlignments;

    // Iterate over the alignments in sorted order
    // If an alignment is distinct (=does not overlap) from the
    // previous alignment, add it to the output collection.
    
    // First alignment is always ok
    outAlignments.push_back(alignments[0]);

    for(size_t i = 1; i < alignments.size(); ++i)
    {
        // Check this alignment against the last alignment added to the output set
        const HapgenAlignment& prevAlign = outAlignments.back();
        const HapgenAlignment& currAlign = alignments[i];

        int s1 = prevAlign.position;
        int e1 = s1 + prevAlign.length;

        int s2 = currAlign.position;
        int e2 = s2 + currAlign.length;
        bool intersecting = Interval::isIntersecting(s1, e1, s2, e2);

        if(prevAlign.referenceID != currAlign.referenceID || !intersecting)
            outAlignments.push_back(currAlign);
    }

    alignments = outAlignments;
}

// Extract a substring for a reference defined by the alignment
std::string HapgenUtil::extractReference(const HapgenAlignment& aln, const ReadTable* pRefTable, int flanking)
{
    const SeqItem& refItem = pRefTable->getRead(aln.referenceID);

    int start = aln.position - flanking;
    int end = aln.position + aln.length + flanking;
    if(start < 0)
        start = 0;
    return refItem.seq.substr(start, end - start);
}


// Extract 3 substrings of a reference described by an alignment
// This writes to 3 out parameters, one for the defined position
// of the reference, one for the upstream flanking bases and one for the downstream flanking
// upstream   defined  downstream
// --------============----------
void HapgenUtil::extractReferenceSubstrings(const HapgenAlignment& aln, 
                                            const ReadTable* pRefTable, 
                                            int flanking,
                                            std::string& outUpstream, 
                                            std::string& outDefined, 
                                            std::string& outDownstream)
{
    const SeqItem& refItem = pRefTable->getRead(aln.referenceID);

    int definedStart = aln.position;
    int definedEnd = aln.position + aln.length;
    int upstreamFlankStart = definedStart - flanking;
    if(upstreamFlankStart < 0)
        upstreamFlankStart = 0;

    if((size_t)definedEnd + flanking > refItem.seq.length())
        flanking = refItem.seq.length() - definedEnd - 1;

    outUpstream = refItem.seq.substr(upstreamFlankStart, definedStart - upstreamFlankStart);
    outDefined = refItem.seq.substr(definedStart, definedEnd - definedStart);
    outDownstream = refItem.seq.substr(definedEnd, flanking);
}

bool HapgenUtil::makeFlankingHaplotypes(const HapgenAlignment& aln, 
                                        const ReadTable* pRefTable, 
                                        int flanking,
                                        const StringVector& inHaplotypes,
                                        StringVector& outHaplotypes)
{
    std::string upstream;
    std::string referenceHaplotype;
    std::string downstream;

    extractReferenceSubstrings(aln, pRefTable, flanking, upstream, referenceHaplotype, downstream);

    // Flip reference strings to match the strand of the input haplotypes
    if(aln.isRC)
    {
        // reverse complement each string
        upstream = reverseComplement(upstream);
        referenceHaplotype = reverseComplement(referenceHaplotype);
        downstream = reverseComplement(downstream);

        // Swap up and downstream
        upstream.swap(downstream);
    }

    // Make the reference haplotype w/ flanking sequence
    std::string referenceFlanking = upstream + referenceHaplotype + downstream;
    outHaplotypes.push_back(referenceFlanking);

    // Check that all sequences match the reference haplotype properly
    bool checkOk = checkAlignmentsAreConsistent(referenceFlanking, inHaplotypes);
    if(!checkOk)
    {
        outHaplotypes.clear();
        return false;
    }

    // Make the flanking sequences for each haplotype
    for(size_t i = 0; i < inHaplotypes.size(); ++i)
    {
        // Skip if the input haplotype exactly matches the reference
        if(inHaplotypes[i] != referenceHaplotype)
            outHaplotypes.push_back(upstream + inHaplotypes[i] + downstream);
    }

    return true;
}
    
// Check that all the strings in the vector align to the same coordinates
// of the passed in sequence
bool HapgenUtil::checkAlignmentsAreConsistent(const std::string& refString, const StringVector& queries)
{
    if(queries.empty())
        return true;

    // Perform local alignments of each query to the refString
    LocalAlignmentResultVector alignments;
    for(size_t i = 0; i < queries.size(); ++i)
        alignments.push_back(StdAlnTools::localAlignment(refString, queries[i]));

    size_t i = 0;
    for(size_t j = 1; j < alignments.size(); ++j)
    {
        if(alignments[i].targetStartPosition != alignments[j].targetStartPosition ||
           alignments[j].targetEndPosition != alignments[j].targetEndPosition)
        {
            std::cerr << "Warning: inconsistent alignments found for haplotype realignment\n";
            std::cerr << "A[" << i << "]: " << alignments[i] << "\n";
            std::cerr << "A[" << j << "]: " << alignments[j] << "\n";
            return false;
        }
    }

    return true;
}

// Extract reads from an FM-index that have a k-mer match to any given haplotypes
// Returns true if the reads were successfully extracted, false if there are 
// more reads than maxReads
bool HapgenUtil::extractHaplotypeReads(const StringVector& haplotypes, 
                                       const BWT* pBWT,
                                       const BWTIntervalCache* pBWTCache,
                                       const SampledSuffixArray* pSSA,
                                       int k,
                                       bool doReverse,
                                       size_t maxReads,
                                       SeqItemVector* pOutReads, 
                                       SeqItemVector* pOutMates)
{
    // Extract the set of reads that have at least one kmer shared with these haplotypes
    // This is a bit of a lengthy procedure with a few steps:
    // 1) extract all the kmers in the haplotypes
    // 2) find the intervals for the kmers in the fm-index
    // 3) compute the set of read indices of the reads from the intervals (using the sampled suffix array)
    // 4) finally, extract the read sequences from the index

    // Make a set of kmers from the haplotypes
    std::set<std::string> kmerSet;
    for(size_t i = 0; i < haplotypes.size(); ++i)
    {
        const std::string& h = haplotypes[i];
        if((int)h.size() < k)
            continue;

        for(size_t j = 0; j < h.size() - k + 1; ++j)
        {
            std::string ks = h.substr(j, k);
            if(doReverse)
                ks = reverseComplement(ks);
            kmerSet.insert(ks);
        }
    }

    // Compute suffix array intervals for the kmers
    size_t totalReads = 0;
    std::vector<BWTInterval> intervals;
    for(std::set<std::string>::const_iterator iter = kmerSet.begin(); iter != kmerSet.end(); ++iter)
    {
        BWTInterval interval = BWTAlgorithms::findIntervalWithCache(pBWT, pBWTCache, *iter);
        intervals.push_back(interval);
        totalReads += interval.size();
    }

    if(totalReads > maxReads)
        return false;

    // Compute the set of reads ids 
    std::set<int64_t> readIndices;
    for(size_t i = 0; i < intervals.size(); ++i)
    {
        BWTInterval interval = intervals[i];
        for(int64_t j = interval.lower; j <= interval.upper; ++j)
        {
            // Get index from sampled suffix array
            SAElem elem = pSSA->calcSA(j, pBWT);
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
        item.seq = BWTAlgorithms::extractString(pBWT, idx);
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
            mateItem.seq = BWTAlgorithms::extractString(pBWT, mateIdx);
            pOutMates->push_back(mateItem);
        }
    }

    return true;
}

// Align a bunch of reads locally to a sequence
LocalAlignmentResultVector HapgenUtil::alignReadsLocally(const std::string& target, const SeqItemVector& reads)
{
    LocalAlignmentResultVector results;
    for(size_t i = 0; i < reads.size(); ++i)
    {
        LocalAlignmentResult fwdAR = StdAlnTools::localAlignment(target, reads[i].seq.toString());
        LocalAlignmentResult rcAR = StdAlnTools::localAlignment(target, reverseComplement(reads[i].seq.toString()));
        results.push_back(fwdAR.score > rcAR.score ? fwdAR : rcAR);
    }
    return results;
}

// Print an alignment to a reference
void HapgenUtil::printAlignment(const std::string& query, const HapgenAlignment& aln, const ReadTable* pRefTable)
{
    const SeqItem& refItem = pRefTable->getRead(aln.referenceID);
    size_t refStart = aln.position;
    size_t refEnd = refStart + aln.length;
    std::string refSubstring = refItem.seq.substr(refStart, refEnd - refStart);

    std::string q = aln.isRC ? reverseComplement(query) : query;
    StdAlnTools::globalAlignment(q, refSubstring, true);
}

