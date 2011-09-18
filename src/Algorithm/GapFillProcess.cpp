///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GapFillProcess - Fill in gaps in a scaffold
//
#include <algorithm>
#include "GapFillProcess.h"
#include "BWTAlgorithms.h"
#include "SGAlgorithms.h"
#include "SGSearch.h"
#include "StdAlnTools.h"
#include "HaplotypeBuilder.h"
#include "MultiAlignment.h"

//
GapFillStats::GapFillStats()
{
    numGapsAttempted = numGapsFilled = 0;
    for(size_t i = 0; i < GFRC_NUM_CODES; ++i)
        numFails[i] = 0;
}

//
void GapFillStats::print() const
{
    printf("Num gaps attempted: %zu\n", numGapsAttempted);
    printf("Num gaps filled: %zu\n", numGapsFilled);
    printf("    Failed -- no anchor: %zu\n", numFails[GFRC_NO_ANCHOR]);
    printf("    Failed -- no haplotype: %zu\n", numFails[GFRC_NO_HAPLOTYPE]);
    printf("    Failed -- ambiguous solution: %zu\n", numFails[GFRC_AMBIGUOUS]);
    printf("    Failed -- assembled size differs from estimate: %zu\n", numFails[GFRC_BAD_SIZE]);
    printf("    Failed -- anchor sequence mismatch: %zu\n", numFails[GFRC_BAD_TRIM]);
}

//
//
//
GapFillProcess::GapFillProcess(const GapFillParameters& params) : m_parameters(params)
{
}

//
GapFillProcess::~GapFillProcess()
{
    m_stats.print();
}

// Process the given scaffold, filling in any gaps found
GapFillResult GapFillProcess::processScaffold(const std::string& scaffold) const
{
    if(m_parameters.verbose > 0)
        std::cout << "Processing scaffold of length " << scaffold.length() << "\n";
    
    GapFillResult result;

    size_t len = scaffold.length();
    size_t currIdx = 0;
    while(currIdx < len)
    {
        char b = scaffold[currIdx];

        if(b != 'N')
        {
            // Append this non-gap character into the output scaffold
            result.scaffold.append(1,b);
            currIdx += 1;
        }
        else
        {
            // Found the start of a gap
            size_t gapLength = 0;
            while(scaffold[currIdx + gapLength] == 'N')
                gapLength += 1;

            if(m_parameters.verbose >= 1)
            printf("Constructing gap at position %zu GapLength: %zu\n", currIdx, gapLength);
            
            // Calculate the left-anchor using the sequence of the scaffold already appended
            AnchorSequence leftAnchor = findAnchor(result.scaffold, result.scaffold.length() - m_parameters.kmer, true);

            // Calculate the right anchor using the sequence of the input scaffold
            AnchorSequence rightAnchor = findAnchor(scaffold, currIdx + gapLength, false);

            // Attempt to build the gap sequence
            std::string gapSequence;
            GapFillReturnCode code = processGap(leftAnchor, rightAnchor, gapSequence);
            m_stats.numGapsAttempted += 1;

            if(code == GFRC_OK)
            {
                // Successfully resolved the gap. Calculate the amount of the anchor sequence
                // that is already present in the scaffold and must be removed. We remove this amount
                // of sequence from the current scaffold, not the gapSequence, to account
                // for the case that the gap sequence is significantly different than the sequence
                // after the left anchor. 
                size_t trimLength = result.scaffold.length() - leftAnchor.position;
                assert(result.scaffold.substr(leftAnchor.position, m_parameters.kmer) == gapSequence.substr(0, m_parameters.kmer));
                result.scaffold.replace(leftAnchor.position, trimLength, gapSequence);
                    
                // We need to update currIdx to point to the next base in the 
                // input scaffold that is not already assembled. This is given
                // by the position of the rightAnchor, plus a kmer
                currIdx = rightAnchor.position + m_parameters.kmer;
                m_stats.numGapsFilled += 1;
            }

            if(code != GFRC_OK)
            {
                // Failed to resolve the gap. Append the gap into the growing scaffold
                m_stats.numFails[code] += 1;

                while(scaffold[currIdx] == 'N')
                {
                    result.scaffold.append(1, 'N');
                    currIdx += 1;
                }
            }
        }
    }

    if(m_parameters.verbose >= 2)
        StdAlnTools::globalAlignment(scaffold, result.scaffold, true);
    return result;
}

// Fill in the specified gap
GapFillReturnCode GapFillProcess::processGap(const AnchorSequence& startAnchor, const AnchorSequence& endAnchor, std::string& outSequence) const
{

    if(m_parameters.verbose > 0)
    {
        std::cout << "\tSTART: " << startAnchor << "\n";
        std::cout << "\tEND: " << endAnchor << "\n";
    }

    if(startAnchor.sequence.empty() || endAnchor.sequence.empty() || startAnchor.sequence == endAnchor.sequence)
        return GFRC_NO_ANCHOR;

    HaplotypeBuilder builder;
    builder.setTerminals(startAnchor, endAnchor);
    builder.setIndex(m_parameters.pBWT, m_parameters.pRevBWT);
    builder.setKmerParameters(m_parameters.kmer, m_parameters.kmerThreshold);
    builder.run();
    
    HaplotypeBuilderResult result;
    builder.parseWalks(result);
    
    if(result.haplotypes.size() > 0)
    {
        if(m_parameters.verbose > 0)
        {
            std::cout << "Built " << result.haplotypes.size() << " candidate haplotypes\n";
            MultiAlignment haplotypeAlignment = MultiAlignmentTools::alignSequences(result.haplotypes);
            haplotypeAlignment.print();
        }

        // Arbitrarily choosing the first haplotype as the sequence
        outSequence = result.haplotypes[0];
        return GFRC_OK;
    }
    else 
    {
        return result.haplotypes.empty() ? GFRC_NO_HAPLOTYPE : GFRC_AMBIGUOUS;
    }
}

AnchorSequence GapFillProcess::findAnchor(const std::string& scaffold, int64_t position, bool upstream) const
{
    WARN_ONCE("TODO: deduplicate findAnchor code");
    AnchorSequence anchor;
    int64_t stride = upstream ? -1 : 1;
    int MAX_DISTANCE = 50;
    int64_t stop = upstream ? position - MAX_DISTANCE : position + MAX_DISTANCE;

    // Cap the travel distance to avoid out of bounds
    if(stop < 0)
        stop = 0;
    if(stop > (int64_t)(scaffold.length() - m_parameters.kmer))
        stop = scaffold.length() - m_parameters.kmer;

    for(; position != stop; position += stride)
    {
        assert(position >= 0);
        std::string testSeq = scaffold.substr(position, m_parameters.kmer);
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
