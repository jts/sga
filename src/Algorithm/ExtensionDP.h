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
// then this alignment can be extended row by row.
#ifndef EXTENSION_DP_H
#define EXTENSION_DP_H

#include <vector>
#include "StdAlnTools.h"

struct DPCell
{
    DPCell() : score(0), ctype('\0') {}
    int score;
    char ctype;
};

typedef std::vector<DPCell> CellVector;

class BandedDPColumn
{
    public:
        BandedDPColumn(int ci, int maxRows, int bandWidth, const BandedDPColumn* prevColumn);

        //
        int getColIdx() const { return m_colIdx; }
        int getRowScore(int row) const;
        char getRowType(int row) const;
        int getCellScore(int colIdx, int rowIdx) const;
        int getMinRow() const { return m_rowStartIdx; }   
        int getMaxRow() const { return m_rowEndIdx; }
        const BandedDPColumn* getPreviousColumn() const;

        //
        void setRowScore(int row, int score, char ctype);
        
        // Calculate the score for the row and set it in the vector
        void fillRowEditDistance(int rowIdx, int matchScore);


    private:

        int getVectorIndex(int rowIdx) const;

        int m_colIdx;
        int m_maxRows;
        int m_rowStartIdx;
        int m_rowEndIdx;
        CellVector m_cells;
        const BandedDPColumn* m_pPrevColumn;
};
typedef std::vector<BandedDPColumn*> BandedDPColumnPtrVector;

namespace ExtensionDP
{
    void initialize(const std::string& extendable, const std::string& fixed, const GlobalAlnParams& params);

    // Calculate the best alignment through the matrix. Assumes that
    // path_t* is pre-allocated with maxPathLength entries.
    void backtrack(const BandedDPColumn* pLastColumn, path_t* path, int* pPathLen, const int maxPathLength);
    
    // Print the alignment starting from the given column
    void printAlignment(const std::string& s1, const std::string& s2, const BandedDPColumn* pLastColumn);

    // Print the full scoring matrix
    void printMatrix(const BandedDPColumnPtrVector& columnPtrVec);
};

#endif
