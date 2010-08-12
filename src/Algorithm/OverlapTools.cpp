//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapTools - Wrapper for the overlap machinery 
// to perform an overlap computation for two strings
//
#include "OverlapTools.h"
#include "SuffixArray.h"
#include "RLBWT.h"

/*
bool OverlapTools::pairOverlap(const std::string& s1, const std::string& s2, 
                               int minOverlap, double errorRate, Overlap& outOverlap)
{
    
    // Construct a read table to hold s1
    ReadTable rt;

    std::string id1 = "id1";
    // Convert the generic string to a DNAString. The DNA string only allows
    // A,C,G,T characters so we permute Ns to one of the bases at random
    DNAString dna1(s1);
    dna1.disambiguate(); 
 
    DNAString dna2(s2);
    dna2.disambiguate();
       
    SeqItem si1 = {id1, dna1};
    rt.addRead(si1);

    // Build the suffix array and RLBWT
    SuffixArray sa(&rt);
    RLBWT bwt(&sa, &rt);



    return false;
}
*/

// Find a suffix of s1 that matches a prefix of s2 with minimal edit distance
void OverlapTools::dpOverlap(const std::string& s1, const std::string& s2)
{
    if(s1.empty() || s2.empty())
        return;

    std::cout << "S1: " << s1.size() << " S2: " << s2.size() << "\n";
    DPAlignment dpAlign(s1, s2, DPM_OVERLAP, DPSS_SIMILARITY);
    //dpAlign.printAlignment(s1, s2);
    findBestOverlapByScore(s1, s2, 1, 0.3, dpAlign);
    /*
    DPMatrix matrix;
    initializeDPMatrixOverlap(m, n, matrix);

    fillDPMatrix(s1, s2, m, n, matrix);
    printDPMatrix(s1, s2, m, n, matrix);
    //printDPAlign(s1, s2, m, n, m - 1, n - 1, matrix);

    // Convert the DP alignment to the maximal overlap match
    findDPBestOverlapSimilarity(s1, s2, m, n, 1, 0.1, matrix);
    */
}    

// Find the cell the best overlap starts in
void OverlapTools::findBestOverlap(const std::string& s1, const std::string& s2, 
                                   int minOverlap, double maxErrorRate, const DPAlignment& dpAlign)
{
    (void)s1;
    (void)s2;
    (void)minOverlap;
    (void)maxErrorRate;
    (void)dpAlign;

    assert(dpAlign.getScoringScheme() == DPSS_EDIT);

    // The last row of the table holds the edit distance for the prefix-suffix overlaps of S1
    // Iterate over this row and calculate the length and edit % of each overlap
    int overlapRowIdx = dpAlign.getNumRows() - 1;
    int numCols = dpAlign.getNumColumns();
    int bestOverlapLen = -1;
    for(int j = minOverlap; j < numCols; ++j)
    {
        int overlapLen = j;

        // Do not allow overlaps to end with insertions or deletions
        // For example, this case is a valid alignment but would be disallowed:
        // S1 AGTACATTTACA--
        // S2    ACATTTACATA
        // This is an overlap of length 11 but the "true" overlap is only 9 bases
        DPOperation lastOp = dpAlign.getPathOperationToCell(s1, s2, overlapRowIdx, j);
        
        if(lastOp != DPO_INSERT)
        {
            int editDistance = dpAlign.getScore(overlapRowIdx,j);
            double errorRate = (double)editDistance / overlapLen;
            if(errorRate < maxErrorRate && overlapLen > bestOverlapLen)
            {
                bestOverlapLen = overlapLen;
            }
            std::cout << "\n";
            dpAlign.printAlignment(s1, s2, overlapRowIdx, j);
            std::cout << "OL: " << overlapLen << " ER: " << errorRate << "\n";
            std::cout << "LO: " << (int)lastOp << " ED: " << editDistance << "\n";
        }
    }
    
    if(bestOverlapLen > 0)
    {
        std::cout << "\nBest overlap: " << bestOverlapLen << "\n";
        dpAlign.printAlignment(s1, s2, overlapRowIdx, bestOverlapLen);
    }
    else
    {
        std::cout << "No significant overlap\n";
    }
}

// Find the cell the best overlap starts in
void OverlapTools::findBestOverlapByScore(const std::string& s1, const std::string& s2, 
                                          int minOverlap, double maxErrorRate, const DPAlignment& dpAlign)
{
    (void)maxErrorRate;
    assert(dpAlign.getScoringScheme() == DPSS_SIMILARITY);

    // The last row of the table holds the edit distance for the prefix-suffix overlaps of S1
    // Iterate over this row and calculate the length and edit % of each overlap
    int overlapRowIdx = dpAlign.getNumRows() - 1;
    int numCols = dpAlign.getNumColumns();
    int bestOverlapLen = -1;
    int bestScore = 0;
    for(int j = minOverlap; j < numCols; ++j)
    {
        int overlapLen = j;

        // Do not allow overlaps to end with insertions or deletions
        // For example, this case is a valid alignment but would be disallowed:
        // S1 AGTACATTTACA--
        // S2    ACATTTACATA
        // This is an overlap of length 11 but the "true" overlap is only 9 bases
        DPOperation lastOp = dpAlign.getPathOperationToCell(s1, s2, overlapRowIdx, j);
        
        //if(lastOp != DPO_INSERT)
        {
            int score = dpAlign.getScore(overlapRowIdx,j);
            if(score > bestScore)
            {
                bestOverlapLen = overlapLen;
                bestScore = score;
            }
            std::cout << "\n";
            dpAlign.printAlignment(s1, s2, overlapRowIdx, j);
            std::cout << "OL: " << overlapLen << "\n";
            std::cout << "LO: " << (int)lastOp << " ED: " << score << "\n";
        }
    }
    
    if(bestOverlapLen > 0)
    {
        std::cout << "\nBest overlap: " << bestOverlapLen << "\n";
        dpAlign.printAlignment(s1, s2, overlapRowIdx, bestOverlapLen);
    }
    else
    {
        std::cout << "No significant overlap\n";
    }
}
