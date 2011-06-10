///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
//
// ExtensionDP - Class implementing an iterative
// dynamic programming alignment. The alignment
// starts by globally constructing an alignment
// between two strings in a seed region, 
// then this alignment can be extended column
// by column for one of the strings.
//
#include "ExtensionDP.h"
#include <stdio.h>
#include <iostream>
#include <assert.h>

static const int MINUS_INF = -1073741823;
inline int max3(int a, int b, int c)
{
    return std::max(std::max(a,b), c);
}

BandedDPColumn::BandedDPColumn(int ci, int maxRows, int bandWidth, const BandedDPColumn* prevColumn)
{
    // Calculate the lower and upper band indices
    m_colIdx = ci;
    m_rowStartIdx = std::max(0, m_colIdx - (bandWidth / 2) - 1);
    m_rowEndIdx = std::min(maxRows, m_colIdx + (bandWidth / 2));

    printf("Col: %d Band: [%d %d]\n", m_colIdx, m_rowStartIdx, m_rowEndIdx);
    int numRows = m_rowEndIdx - m_rowStartIdx + 1;
    m_scores.resize(numRows, MINUS_INF);

    m_pPrevColumn = prevColumn;
}

// Get the score for row
int BandedDPColumn::getRowScore(int row) const
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx == -1)
        return MINUS_INF;
    else
        return m_scores[vec_idx];
}

// Set the score for row
void BandedDPColumn::setRowScore(int row, int score)
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx != -1)
        m_scores[vec_idx] = score;
}

int BandedDPColumn::getCellScore(int colIdx, int rowIdx) const
{
    // If accessing the col/row with negative index, return
    // the initial gap penalties. For now, we allow free gaps
    // for the first row
    if(colIdx == -1 || rowIdx == -1)
        return 0;

    const BandedDPColumn* pColumn;
    if(colIdx < m_colIdx)
        pColumn = m_pPrevColumn;
    else
        pColumn = this;
    assert(pColumn != NULL);
    
    return pColumn->getRowScore(rowIdx);   
}

// Calculate the score for the given row in the column
void BandedDPColumn::fillRow(int rowIdx, int matchScore, const GlobalAlnParams& params)
{
    int vec_idx = getVectorIndex(rowIdx);
    if(vec_idx == -1)
        return; // out of band, ignore

    // Get the scores for the cell above, to the left and diagonal
    int above = getCellScore(m_colIdx, rowIdx - 1) - params.gap_open;
    int diag = getCellScore(m_colIdx - 1, rowIdx - 1) + matchScore;
    int left = getCellScore(m_colIdx - 1, rowIdx) - params.gap_open;
    int score = max3(above, left, diag);
    printf("Cell(%d, %d) a: %d d: %d l: %d s: %d\n", m_colIdx, rowIdx, above, diag, left, score);
    setRowScore(rowIdx, score);
}

// Transform a row index into an index in the scores vector
int BandedDPColumn::getVectorIndex(int rowIdx) const
{
    if(rowIdx < m_rowStartIdx || rowIdx > m_rowEndIdx)
        return -1;
    else
        return rowIdx - m_rowStartIdx;
}

// Initialize the extension DP by computing a global alignment between extendable and fixed
void ExtensionDP::initialize(const std::string& extendable, const std::string& fixed, const GlobalAlnParams& params)
{
    std::vector<BandedDPColumn*> m_columnPtrs;
    
    // Initialize the zero-th column
    size_t numRows = fixed.size() + 1;
    size_t numCols = extendable.size() + 1;
    std::cout << "ROWS: " << numRows << "\n";
    BandedDPColumn* pZeroCol = new BandedDPColumn(0, numRows, params.bandwidth, NULL);

    // Set the score of (0,0) to be zero
    pZeroCol->setRowScore(0, 0);
    for(size_t i = 1; i < numRows; ++i)
        pZeroCol->setRowScore(i, -(params.gap_open + i*params.gap_ext));

    m_columnPtrs.push_back(pZeroCol);

    for(size_t colIdx = 1; colIdx < numCols; ++colIdx)
    {
        BandedDPColumn* pPrevCol = m_columnPtrs[colIdx - 1];
        BandedDPColumn* pCurrCol = new BandedDPColumn(colIdx, numRows, params.bandwidth, pPrevCol);

        // Set the first row score
        pCurrCol->setRowScore(0, -(params.gap_open + colIdx*params.gap_ext));
        
        // Fill in the rest of the row
        int startRow = std::max(1, pCurrCol->getMinRow());
        for(int rowIdx = startRow; rowIdx <= pCurrCol->getMaxRow(); ++rowIdx)
        {
            int match_score = extendable[colIdx - 1] == fixed[rowIdx - 1] ? params.match : -params.mismatch;
            pCurrCol->fillRow(rowIdx, match_score, params);
        }
        m_columnPtrs.push_back(pCurrCol);
    }

    printMatrix(m_columnPtrs);

    std::cout << "STDALN: \n";
    StdAlnTools::printGlobalAlignment(extendable, fixed);

    // Delete the columns
    for(size_t i = 0; i < m_columnPtrs.size(); ++i)
        delete m_columnPtrs[i];
}

void ExtensionDP::printMatrix(const BandedDPColumnPtrVector& columnPtrVec)
{
    int numCols = columnPtrVec.size();
    int numRows = columnPtrVec[numCols - 1]->getMaxRow();

    for(int i = 0; i < numRows; ++i)
    {
        for(int j = 0; j < numCols; ++j)
        {
            int score = columnPtrVec[j]->getRowScore(i);
            printf("%d\t", score);
        }
        printf("\n");
    }
}

