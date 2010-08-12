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

    size_t m = s1.size() + 1; // = num rows
    size_t n = s2.size() + 1; // = num columns

    DPMatrix matrix;
    initializeDPMatrixOverlap(m, n, matrix);

    fillDPMatrix(s1, s2, m, n, matrix);
    printDPMatrix(s1, s2, m, n, matrix);
    //printDPAlign(s1, s2, m, n, m - 1, n - 1, matrix);

    // Convert the DP alignment to the maximal overlap match
    findDPBestOverlap(s1, s2, m, n, 1, 0.1, matrix); 
}    

// Setup the dynamic programming matrix for overlap computation
void OverlapTools::initializeDPMatrixOverlap(int m, int n, DPMatrix& score)
{
    score.resize(m); // rows
    for(int i = 0; i < m; ++i)
        score[i].resize(n); // columns

    // Set score[i][0] (first column) to be zero for all i
    // This is so we don't penalize the portion of S1 that
    // S2 does not overlap
    for(int i = 0; i < m; ++i)
        score[i][0] = 0;

    // Set score [0][j] (first row) to be j for all j
    // This represents not matching the first j bases of S1
    for(int j = 0; j < n; ++j)
        score[0][j] = j;
}

// Fill in the DP matrix
void OverlapTools::fillDPMatrix(const std::string& s1, const std::string& s2, int m, int n, DPMatrix& score)
{
    // Update the value of score[i][j] for all i and j greater than zero
    for(int i = 1; i < m; ++i)
    {
        for(int j = 1; j < n; ++j)
        {
            int up = score[i - 1][j] + 1; // skip a base of S1
            int left = score[i][j - 1] + 1; // skip a base of S2

            int idx1 = i - 1;
            int idx2 = j - 1;
            int match_score = (s1[idx1] == s2[idx2] ? 0 : 1);
            int diag = score[i - 1][j - 1] + match_score;

            score[i][j] = min3(up, left, diag);
        }
    }
}

// Find the cell the best overlap starts in
void OverlapTools::findDPBestOverlap(const std::string& s1, const std::string& s2, 
                                     int m, int n, int minOverlap, double maxErrorRate,
                                     const DPMatrix& matrix)
{
    // The last row of the table holds the edit distance for the prefix-suffix overlaps of S1
    // Iterate over this row and calculate the length and edit % of each overlap
    int overlapRowIdx = m - 1;
    int bestOverlapLen = -1;
    int prevEditDistance = (minOverlap > 0) ? matrix[overlapRowIdx][minOverlap - 1] : 0;
    for(int j = minOverlap; j < n; ++j)
    {
        int overlapLen = j;
        int editDistance = matrix[overlapRowIdx][j];
        double errorRate = (double)editDistance / overlapLen;

        // Do not allow overlaps to end with insertions or deletions
        // For example, this case would be disallowed:
        // S1 AGTACATTTACA--
        // S2    ACATTTACATA
        // This is an overlap of length 11 but the "true" overlap is only 9 bases
        bool hasTerminalInsertion = (editDistance == prevEditDistance + 1);
        //if(!hasTerminalInsertion)
        {
            if(errorRate < maxErrorRate && overlapLen > bestOverlapLen)
            {
                bestOverlapLen = overlapLen;
            }
            std::cout << "\n";
            printDPAlign(s1, s2, m, n, overlapRowIdx, j, matrix);
            std::cout << "OL: " << overlapLen << " ER: " << errorRate << "\n";
            std::cout << "HT: " << hasTerminalInsertion << " ED: " << editDistance << " prev: " << prevEditDistance << "\n";
        }
        prevEditDistance = editDistance;
    }
    
    if(bestOverlapLen > 0)
    {
        std::cout << "\nBest overlap: " << bestOverlapLen << "\n";
        printDPAlign(s1, s2, m, n, overlapRowIdx, bestOverlapLen, matrix);
    }
    else
    {
        std::cout << "No significant overlap\n";
    }
}

OverlapTools::DPPath OverlapTools::calculateDPPath(const std::string& s1, const std::string& s2, 
                         int m, int n, int startI, int startJ, 
                         const DPMatrix& scores)
{
    DPPath out;
    DPPathNode currNode;
    currNode.i = startI;
    currNode.j = startJ;
    currNode.op = s1[m - 2] == s2[n - 2] ? DPO_MATCH : DPO_SUB;

    while(currNode.i != 0 || currNode.j != 0)
    {
        // Choose the cell that preceded the current cell
        // in the optimal path. If multiple cells are valid,
        // choose the diagonal
        DPPathNode nextNode;
        bool found = false;
        int currScore = scores[currNode.i][currNode.j];
        //printf("currNode: (%d %d)\n", currNode.i, currNode.j);

        // Check diagonal first
        int i = currNode.i - 1;
        int j = currNode.j - 1;
        if(i >= 0 && j >= 0)
        {
            int match_score = (s1[i] == s2[j] ? 0 : 1);
            if(currScore == scores[i][j] + match_score)
            {
                nextNode.i = i;
                nextNode.j = j;
                nextNode.op = (s1[i] == s2[j]) ? DPO_MATCH : DPO_SUB;
                found = true;
            }
        }

        // Check above element
        i = currNode.i - 1;
        j = currNode.j;
        if(!found && i >= 0 && j >= 0)
        {
            // In the first column, there is no penalty
            // for a deletion
            int above_score = (j == 0) ? 0 : 1;
            if(currScore == scores[i][j] + above_score)
            {
                nextNode.i = i;
                nextNode.j = j;
                nextNode.op = DPO_DELETE;
                found = true;
            }
        }

        // Check left element
        i = currNode.i;
        j = currNode.j - 1;
        if(!found && i >= 0 && j >= 0)
        {                  
            if(currScore == scores[i][j] + 1)
            {
                nextNode.i = i;
                nextNode.j = j;
                nextNode.op = DPO_INSERT;
                found = true;
            }
        }

        assert(found);
        out.push_back(nextNode);
        currNode = nextNode;
    }
    out.reverse();
    return out;
}


void OverlapTools::printDPMatrix(const std::string& s1, const std::string& s2, int m, int n, const DPMatrix& score)
{
    // Print the DP matrix
    std::cout << " \t-\t";
    for(int j = 1; j < n; ++j)
    {
        std::cout << s2[j - 1] << "\t";
    }
    std::cout << "\n";

    for(int i = 0; i < m; ++i)
    {
        if(i == 0)
            std::cout << "-\t";
        else
            std::cout << s1[i - 1] << "\t";

        for(int j = 0; j < n; ++j)
        {
            std::cout << score[i][j] << "\t";
        }
        std::cout << "\n";
    }

    printDPAlign(s1, s2, m, n, m - 1, n - 1, score);
}

// Print the best alignment indicated by the matrix starting from (startI, startJ)
void OverlapTools::printDPAlign(const std::string& s1, const std::string& s2, 
                                int m, int n, int startI, int startJ, 
                                const DPMatrix& scores)
{
    if(scores.empty())
        return;
    
    int edit_distance = scores[startI][startJ];
    DPPath backtrack = calculateDPPath(s1, s2, m, n, startI, startJ, scores);

    // Print the match strings
    std::string m1;
    std::string m2;
    std::string os;
    for(DPPath::iterator iter = backtrack.begin(); iter != backtrack.end(); ++iter)
    {
        int i = iter->i;
        int j = iter->j;
        //printf("(%d %d)\n", i, j);
        switch(iter->op)
        {
            case DPO_MATCH:
            case DPO_SUB:
                assert(i >= 0 && i < (int)s1.size() && j >= 0 && j < (int)s2.size());
                m1.push_back(s1[i]);
                m2.push_back(s2[j]);
                os.push_back(iter->op == DPO_MATCH ? '.' : 'S');
                break;
            case DPO_INSERT:
                assert(j >= 0 && j < (int)s2.size());
                m1.push_back('-');
                m2.push_back(s2[j]);
                os.push_back('I');
                break;
            case DPO_DELETE:
                assert(i >= 0 && i < (int)s1.size());
                m1.push_back(s1[i]);
                m2.push_back('-');
                os.push_back('D');
                break;
        }
    }

    std::cout << "M1: " << m1 << "\n";
    std::cout << "M2: " << m2 << "\n";
    std::cout << "OS: " << os << "\n";
    std::cout << "Edit distance: " << edit_distance << "\n";
}
