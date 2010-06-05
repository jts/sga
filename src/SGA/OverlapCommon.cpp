//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapCommon - Common wrapper used for finding overlaps
// for a set of reads
//
#include "OverlapCommon.h"

// Convert a line from a hits file into a vector of overlaps and sets the flag
// indicating whether the read was found to be a substring of other reads
// Only the forward read table is used since we only care about the IDs and length
// of the read, not the sequence, so that we don't need an explicit reverse read table
void OverlapCommon::parseHitsString(const std::string& hitString, 
                                    const ReadInfoTable* pRIT, 
                                    const SuffixArray* pFwdSAI, const SuffixArray* pRevSAI, 
                                    size_t& readIdx, OverlapVector& outVector, bool& isSubstring)
{
    OverlapVector outvec;
    std::istringstream convertor(hitString);

    // Read the overlap blocks for a read
    size_t numBlocks;
    convertor >> readIdx >> isSubstring >> numBlocks;

    //std::cout << "<Read> idx: " << readIdx << " count: " << numBlocks << "\n";
    for(size_t i = 0; i < numBlocks; ++i)
    {
        // Read the block
        OverlapBlock record;
        convertor >> record;
        //std::cout << "\t" << record << "\n";

        // Iterate through the range and write the overlaps
        for(int64_t j = record.ranges.interval[0].lower; j <= record.ranges.interval[0].upper; ++j)
        {
            const SuffixArray* pCurrSAI = (record.flags.isTargetRev()) ? pRevSAI : pFwdSAI;
            const ReadInfo& queryInfo = pRIT->getReadInfo(readIdx);

            int64_t saIdx = j;

            // The index of the second read is given as the position in the SuffixArray index
            const ReadInfo& targetInfo = pRIT->getReadInfo(pCurrSAI->get(saIdx).getID());

            // Skip self alignments and non-canonical (where the query read has a lexo. higher name)
            if(queryInfo.id != targetInfo.id)
            {    
                // Compute the endpoints of the overlap
                int s1 = queryInfo.length - record.overlapLen;
                int e1 = s1 + record.overlapLen - 1;
                SeqCoord sc1(s1, e1, queryInfo.length);

                int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
                int e2 = s2 + record.overlapLen - 1;
                SeqCoord sc2(s2, e2, targetInfo.length);

                // The coordinates are always with respect to the read, so flip them if
                // we aligned to/from the reverse of the read
                if(record.flags.isQueryRev())
                    sc1.flip();
                if(record.flags.isTargetRev())
                    sc2.flip();

                bool isRC = record.flags.isTargetRev() != record.flags.isQueryRev();

                Overlap o(queryInfo.id, sc1, targetInfo.id, sc2, isRC, record.numDiff);
            
                // The alignment logic above has the potential to produce duplicate alignments
                // To avoid this, we skip overlaps where the id of the first coord is lexo. lower than 
                // the second or the match is a containment and the query is reversed (containments can be 
                // output up to 4 times total).
                if(o.id[0] < o.id[1] || (o.match.isContainment() && record.flags.isQueryRev()))
                    continue;

                outVector.push_back(o);
            }
        }
    }
}
