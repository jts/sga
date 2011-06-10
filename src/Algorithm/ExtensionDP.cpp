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
static const int POSITIVE_INF = 1073741823;

inline int max3(int a, int b, int c)
{
    return std::max(std::max(a,b), c);
}

inline int min3(int a, int b, int c)
{
    return std::min(std::min(a,b), c);
}

BandedDPColumn::BandedDPColumn(int ci, int maxRows, int bandWidth, const BandedDPColumn* prevColumn)
{
    // Calculate the lower and upper band indices
    m_colIdx = ci;
    m_rowStartIdx = std::max(0, m_colIdx - (bandWidth / 2) - 1);
    m_rowEndIdx = std::min(maxRows - 1, m_colIdx + (bandWidth / 2));

    int numRows = m_rowEndIdx - m_rowStartIdx + 1;
    m_cells.resize(numRows);

    m_pPrevColumn = prevColumn;
}

// Get the score for row
int BandedDPColumn::getRowScore(int row) const
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx == -1)
        return POSITIVE_INF;
    else
        return m_cells[vec_idx].score;
}

// Get the type of the row
char BandedDPColumn::getRowType(int row) const
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx == -1)
        return FROM_M;
    else
        return m_cells[vec_idx].ctype;
}

//
int BandedDPColumn::getBestRowIndex() const
{
    int bestScore = POSITIVE_INF;
    int bestIdx = -1;
    for(int i = m_rowStartIdx; i <= m_rowEndIdx; ++i)
    {
        int s = getRowScore(i);
        if(s < bestScore)
        {
            bestScore = s;
            bestIdx = i;
        }
    }
    return bestIdx;
}

// Return the pointer to the previous column
const BandedDPColumn* BandedDPColumn::getPreviousColumn() const
{
    return m_pPrevColumn;
}

// Set the score for row
void BandedDPColumn::setRowScore(int row, int score, char ctype)
{
    int vec_idx = getVectorIndex(row);
    if(vec_idx != -1)
    {
        m_cells[vec_idx].score = score;
        m_cells[vec_idx].ctype = ctype;
    }
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

// Calculate the score for the given row in the column using an edit distance scheme
void BandedDPColumn::fillRowEditDistance(int rowIdx, int matchScore)
{
    int vec_idx = getVectorIndex(rowIdx);
    if(vec_idx == -1)
        return; // out of band, ignore

    // Get the scores for the cell above, to the left and diagonal
    int above = getCellScore(m_colIdx, rowIdx - 1) + 1;
    int diag = getCellScore(m_colIdx - 1, rowIdx - 1) + matchScore;
    int left = getCellScore(m_colIdx - 1, rowIdx) + 1;
    int score = min3(above, left, diag);

    char ctype;
    if(score == above)
        ctype = FROM_D;
    else if(score == diag)
        ctype = FROM_M;
    else
        ctype = FROM_I;
//    printf("Cell(%d, %d) a: %d d: %d l: %d s: %d t: %d\n", m_colIdx, rowIdx, above, diag, left, score, ctype);
    setRowScore(rowIdx, score, ctype);
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
// This function allocates memory and stores the created pointers in outPtrVec, which
// should initially be empty
void ExtensionDP::createInitialAlignment(const std::string& extendable, const std::string& fixed, int bandwidth, BandedDPColumnPtrVector& outPtrVec)
{
    assert(outPtrVec.empty());
    // Initialize the zero-th column
    size_t numRows = fixed.size() + 1;
    size_t numCols = extendable.size() + 1;

    // Set the score of column zero 
    BandedDPColumn* pZeroCol = new BandedDPColumn(0, numRows, bandwidth, NULL);
    pZeroCol->setRowScore(0, 0, FROM_M);
    for(size_t i = 1; i < numRows; ++i)
        pZeroCol->setRowScore(i, i, FROM_I);

    outPtrVec.push_back(pZeroCol);

    for(size_t colIdx = 1; colIdx < numCols; ++colIdx)
    {
        BandedDPColumn* pPrevCol = outPtrVec[colIdx - 1];
        BandedDPColumn* pCurrCol = new BandedDPColumn(colIdx, numRows, bandwidth, pPrevCol);

        // Set the first row score
        pCurrCol->setRowScore(0, colIdx, FROM_D);
        
        // Fill in the rest of the row
        int startRow = std::max(1, pCurrCol->getMinRow());
        for(int rowIdx = startRow; rowIdx <= pCurrCol->getMaxRow(); ++rowIdx)
        {
            int match_score = extendable[colIdx - 1] == fixed[rowIdx - 1] ? 0 : 1;
            pCurrCol->fillRowEditDistance(rowIdx, match_score);
        }
        outPtrVec.push_back(pCurrCol);
    }
}

//
//
double ExtensionDP::calculateLocalEditPercentage(const BandedDPColumn* pColumn, int numBases)
{
    // Get the row with in the start column with the best score
    int rowIdx = pColumn->getBestRowIndex();
    assert(rowIdx != -1);
    int colIdx = pColumn->getColIdx();
    int startScore = pColumn->getRowScore(rowIdx);

    printf("Start cell(%d %d) = %d\n", rowIdx, colIdx, startScore);
    int alignLength = 1;
    while(alignLength < numBases)
    {
        char ctype = pColumn->getRowType(rowIdx);
        switch(ctype)
        {
   			case FROM_M: 
                pColumn = pColumn->getPreviousColumn();
                --colIdx; 
                --rowIdx; 
                break;
			case FROM_D: 
                --rowIdx; 
                break;
			case FROM_I: 
                pColumn = pColumn->getPreviousColumn();
                --colIdx; 
                break;
        }
        alignLength += 1;
    }
    int endScore = pColumn->getRowScore(rowIdx);
    printf("End cell(%d %d) = %d al: %d\n", rowIdx, colIdx, endScore, alignLength);

    assert(startScore >= endScore);
    int numDiffs = startScore - endScore;
    return static_cast<double>(numDiffs) / alignLength;
}

// Calculate the best alignment through the matrix. Assumes that
// path_t* is pre-allocated with max_path entries.
void ExtensionDP::backtrack(const BandedDPColumn* pLastColumn, path_t* path, int* path_len, const int maxPathLength)
{
    const BandedDPColumn* pCurrColumn = pLastColumn;

    // Calculate the starting row and column
    int colIdx = pLastColumn->getColIdx();
    int rowIdx = pLastColumn->getMaxRow();
    int pathLength = 0;


    path_t* p = &path[pathLength];
    char ctype = pLastColumn->getRowType(rowIdx);

	p->ctype = ctype; 
    p->i = rowIdx; 
    p->j = colIdx;
	pathLength += 1;
    
	do {
        assert(pathLength < maxPathLength);
        assert(pCurrColumn != NULL);

        p = &path[pathLength];

        // Update indices
		switch (ctype) {
			case FROM_M: 
                pCurrColumn = pCurrColumn->getPreviousColumn();
                --colIdx; 
                --rowIdx; 
                break;
			case FROM_D: 
                --rowIdx; 
                break;
			case FROM_I: 
                pCurrColumn = pCurrColumn->getPreviousColumn();
                --colIdx; 
                break;
		}
		ctype = pCurrColumn->getRowType(rowIdx);
		p->ctype = ctype; 
        p->i = rowIdx; 
        p->j = colIdx;
        pathLength += 1;
	} while (rowIdx || colIdx);

    // We ignore the path character from the first column/row, like stdaln
    *path_len = pathLength - 1;
}

void ExtensionDP::printAlignment(const std::string& extendable, const std::string& fixed, const BandedDPColumn* pLastColumn)
{
    // allocate space for the path
    int maxPathLength = extendable.size() + fixed.size();
    path_t* path = new path_t[maxPathLength];
    int path_len;

    // Perform the backtrack
    backtrack(pLastColumn, path, &path_len, maxPathLength);

    // Print the strings
    std::string paddedF, paddedE, paddedMatch;
    StdAlnTools::makePaddedStrings(fixed, extendable, path, path_len, paddedF, paddedE, paddedMatch);
    StdAlnTools::printPaddedStrings(paddedF, paddedE, paddedMatch);

    std::cout << "CIGAR: " << StdAlnTools::makeCigar(path, path_len) << "\n";

    delete [] path;
}

void ExtensionDP::printMatrix(const BandedDPColumnPtrVector& columnPtrVec)
{
    int numCols = columnPtrVec.size();
    int numRows = columnPtrVec[numCols - 1]->getMaxRow() + 1;

    printf("c:\t");
    for(int j = 0; j < numCols; ++j)
    {
        printf("%d\t", j);
    }
    printf("\n");

    for(int i = 0; i < numRows; ++i)
    {
        printf("%d:\t", i);
        for(int j = 0; j < numCols; ++j)
        {
            int score = columnPtrVec[j]->getRowScore(i);
            printf("%d\t", score);
        }
        printf("\n");
    }
}

