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

    // Find gaps in the scaffold and attempt to build sequence for them
    int gapStart = -1;
    int gapEnd = -1;
    for(size_t i = 0; i < scaffold.length(); ++i)
    {
        char b = scaffold[i];
        if(b == 'N')
        {
            if(gapStart == -1)
                gapStart = i; // start a new gap
            gapEnd = i;
        }
        else
        {
            // end a gap
            if(gapStart != -1)
            {
                std::string gapSequence;
                GapFillReturnCode code = processGap(scaffold, gapStart, gapEnd, gapSequence);
                m_stats.numGapsAttempted += 1;

                if(code == GFRC_OK)
                {
                    result.scaffold.append(gapSequence);
                    m_stats.numGapsFilled += 1;
                }
                else
                {
                    m_stats.numFails[code] += 1;
                    result.scaffold.append(gapEnd - gapStart + 1, 'N');
                }
                
                // Reset gap positions
                gapStart = gapEnd = -1;
            }

            // Append in the normal character to the scaffold
            result.scaffold.append(1, b);
        }
    }

    //StdAlnTools::globalAlignment(scaffold, outScaffold, true);
    return result;
}

// Fill in the specified gap
GapFillReturnCode GapFillProcess::processGap(const std::string& scaffold, int gapStart, int gapEnd, std::string& outSequence) const
{
    if(m_parameters.verbose > 0)
        std::cout << "Attempting to fill gap of length " << gapEnd - gapStart + 1 << "\n";
    AnchorSequence startAnchor = findAnchor(scaffold, gapStart - m_parameters.kmer, true);
    AnchorSequence endAnchor = findAnchor(scaffold, gapEnd + 1, false);

    if(startAnchor.sequence.empty() || endAnchor.sequence.empty() || startAnchor.sequence == endAnchor.sequence)
        return GFRC_NO_ANCHOR;

    if(m_parameters.verbose > 0)
    {
        std::cout << "\tSTART: " << startAnchor << "\n";
        std::cout << "\tEND: " << endAnchor << "\n";
    }

    HaplotypeBuilder builder;
    builder.setTerminals(startAnchor, endAnchor);
    builder.setIndex(m_parameters.pBWT, m_parameters.pRevBWT);
    builder.setKmerParameters(m_parameters.kmer, m_parameters.kmerThreshold);
    builder.run();
    
    HaplotypeBuilderResult result;
    builder.parseWalks(result);
    
    
    if(result.haplotypes.size() >= 1)
    {
        if(m_parameters.verbose > 0)
        {
            std::cout << "Built " << result.haplotypes.size() << " candidate haplotypes\n";
            MultiAlignment haplotypeAlignment = MultiAlignmentTools::alignSequences(result.haplotypes);
            haplotypeAlignment.print();
        }

        // Arbitrarily choosing the first haplotype as the sequence
        std::string& selectedSeq = result.haplotypes[0];
        
        // Map the coordinates of the gap start/end onto the haplotype
        assert(gapStart >= startAnchor.position);
        assert(gapEnd <= endAnchor.position);
        int hap_gap_start = (gapStart - startAnchor.position);
        int hap_gap_end  = selectedSeq.length() - m_parameters.kmer - (endAnchor.position - gapEnd);

        //printf("GapStart: %d GapEnd: %d AnchorStart: %d AnchorEnd: %d\n", gapStart, gapEnd, startAnchor.position, endAnchor.position);
        //std::cout << "GSL: " << selectedSeq.length() << " HS: " << hap_gap_start << " HE: " << hap_gap_end << "\n";

        // Extract the sequence of the gap
        outSequence = selectedSeq.substr(hap_gap_start, hap_gap_end - hap_gap_start + 1);
    }
    
    return !result.haplotypes.empty() ? GFRC_OK : GFRC_NO_HAPLOTYPE;
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
