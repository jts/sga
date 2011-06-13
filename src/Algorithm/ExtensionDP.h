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
        BandedDPColumn(int ci, int maxRows, int bandwidth, const BandedDPColumn* prevColumn);

        //
        int getColIdx() const { return m_colIdx; }
        int getBandwidth() const { return m_bandwidth; }

        int getRowScore(int row) const;
        char getRowType(int row) const;
        int getCellScore(int colIdx, int rowIdx) const;
        int getMinRow() const { return m_rowStartIdx; }   
        int getMaxRow() const { return m_rowEndIdx; }
        const BandedDPColumn* getPreviousColumn() const;

        //
        void setRowScore(int row, int score, char ctype);
        
        // Returns the row index of the best scoring cell
        int getBestRowIndex() const;

        // Calculate the score for the row and set it in the vector
        void fillRowEditDistance(int rowIdx, int matchScore);

        // Returns true if the best alignment in this column is to the endpoint of the query sequence
        bool isAlignedToEnd() const;

    private:

        int getVectorIndex(int rowIdx) const;

        int m_colIdx;
        int m_maxRows;
        int m_rowStartIdx;
        int m_rowEndIdx;
        int m_bandwidth;
        CellVector m_cells;
        const BandedDPColumn* m_pPrevColumn;
};
typedef std::vector<BandedDPColumn*> BandedDPColumnPtrVector;

namespace ExtensionDP
{
    // Initialize the extension DP by computing a global alignment between extendable and fixed
    // This function allocates memory and stores the created pointers in outPtrVec
    void createInitialAlignment(const std::string& extendable, const std::string& fixed, int bandwidth, BandedDPColumnPtrVector& outPtrVec);

    // Create a new column, representing the extended alignment from pPrevColumn to char b
    BandedDPColumn* createNewColumn(char b, const std::string& fixed, const BandedDPColumn* pPrevColumn);

    // Calculate the best alignment through the matrix. Assumes that
    // path_t* is pre-allocated with maxPathLength entries.
    void backtrack(const BandedDPColumn* pLastColumn, path_t* path, int* pPathLen, const int maxPathLength);

    // Calculate the edit percentage of the alignment starting from the given column, over the last numBases
    double calculateLocalEditPercentage(const BandedDPColumn* pStartColumn, int numBases);
    double calculateGlobalEditPercentage(const BandedDPColumn* pStartColumn);
    void countEditsAndAlignLength(const BandedDPColumn* pStartColumn, int& edits, int& alignLength);

    // Print the alignment starting from the given column
    void printAlignment(const std::string& s1, const std::string& s2, const BandedDPColumn* pLastColumn);

    // Print the full scoring matrix
    void printMatrix(const BandedDPColumnPtrVector& columnPtrVec);
};

#endif
