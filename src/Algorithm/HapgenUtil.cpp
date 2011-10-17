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

// Align the haplotype to the reference genome represented by the BWT/SSA pair
HapgenAlignmentVector HapgenUtil::alignHaplotypeToReference(const std::string& haplotype,
                                                            const BWT* pReferenceBWT,
                                                            const SampledSuffixArray* pReferenceSSA)
{
    HapgenAlignmentVector results;

    LRAlignment::LRParams params;
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
                HapgenAlignment aln(hits[j].targetID, hits[j].t_start, hits[j].length, i == 1);
                results.push_back(aln);
            }
            else
            {
                std::cerr << "Skipped partial alignment (" << q_alignment_length << " of " << haplotype.length() << ")\n";
            }
        }
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

