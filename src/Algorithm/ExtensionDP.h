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

typedef std::vector<int> IntVector;

class BandedDPColumn
{
    public:
        BandedDPColumn(int ci, int maxRows, int bandWidth, const BandedDPColumn* prevColumn);

        //
        int getRowScore(int row) const;
        int getCellScore(int colIdx, int rowIdx) const;
        int getMinRow() const { return m_rowStartIdx; }   
        int getMaxRow() const { return m_rowEndIdx; }

        //
        void setRowScore(int row, int score);
        
        // Calculate the score for the row and set it in the vector
        void fillRowEditDistance(int rowIdx, int matchScore);

    private:

        int getVectorIndex(int rowIdx) const;

        int m_colIdx;
        int m_maxRows;
        int m_rowStartIdx;
        int m_rowEndIdx;
        IntVector m_scores;
        const BandedDPColumn* m_pPrevColumn;
};
typedef std::vector<BandedDPColumn*> BandedDPColumnPtrVector;

namespace ExtensionDP
{
    void initialize(const std::string& extendable, const std::string& fixed, const GlobalAlnParams& params);
    void printMatrix(const BandedDPColumnPtrVector& columnPtrVec);

};

#endif
