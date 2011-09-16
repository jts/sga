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
//
//
GapFillProcess::GapFillProcess(const GapFillParameters& params) : m_parameters(params)
{
    m_gapsAttempted = 0;
    m_gapsFilled = 0;
}

//
GapFillProcess::~GapFillProcess()
{
    std::cout << "Gaps attempted: " << m_gapsAttempted << "\n";
    std::cout << "Gaps filled: " << m_gapsFilled << "\n";
}

// Process the given scaffold, filling in any gaps found
void GapFillProcess::processScaffold(const std::string& scaffold) const
{
    if(m_parameters.verbose > 0)
        std::cout << "Processing scaffold of length " << scaffold.length() << "\n";
    
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
                bool success = processGap(scaffold, gapStart, gapEnd);
                m_gapsAttempted += 1;
                if(success)
                    m_gapsFilled += 1;

                gapStart = gapEnd = -1;
            }
        }
    }
}

// Fill in the specified gap
bool GapFillProcess::processGap(const std::string& scaffold, int gapStart, int gapEnd) const
{
    if(m_parameters.verbose > 0)
        std::cout << "Attempting to fill gap of length " << gapEnd - gapStart + 1 << "\n";
    AnchorSequence startAnchor = findAnchor(scaffold, gapStart - m_parameters.kmer, true);
    AnchorSequence endAnchor = findAnchor(scaffold, gapEnd, false);
    if(startAnchor.sequence.empty() || endAnchor.sequence.empty() || startAnchor.sequence == endAnchor.sequence)
        return false;

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
    
    
    if(result.haplotypes.size() >= 1 && m_parameters.verbose > 0)
    {
        std::cout << "Built " << result.haplotypes.size() << " candidate haplotypes\n";
        MultiAlignment haplotypeAlignment = MultiAlignmentTools::alignSequences(result.haplotypes);
        haplotypeAlignment.print();
    }
    
    return !result.haplotypes.empty();
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
