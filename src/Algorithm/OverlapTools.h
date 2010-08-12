//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapTools - Wrapper for the overlap machinery 
// to perform an overlap computation for two strings
//
#ifndef OVERLAPTOOLS_H
#define OVERLAPTOOLS_H

#include "Match.h"
#include <list>

namespace OverlapTools
{
    // Datatypes
    typedef std::vector<int> IntVector;
    typedef std::vector<IntVector> DPMatrix;

    enum DPOperation
    {
        DPO_MATCH,
        DPO_SUB,
        DPO_INSERT,
        DPO_DELETE
    };

    struct DPPathNode
    {
        DPPathNode() : i(0), j(0), op(DPO_MATCH) {}
        int i;
        int j;
        DPOperation op;
    };

    typedef std::list<DPPathNode> DPPath;

    inline int min3(int a, int b, int c)
    {
        return std::min(std::min(a,b), c);
    }

    void dpOverlap(const std::string& s1, const std::string& s2);
    void initializeDPMatrixOverlap(int m, int n, DPMatrix& score);
    void fillDPMatrix(const std::string& s1, const std::string& s2, int m, int n, DPMatrix& score);
    void findDPBestOverlap(const std::string& s1, const std::string& s2, 
                                     int m, int n, int minOverlap, double maxErrorRate,
                                     const DPMatrix& matrix);

    DPPath calculateDPPath(const std::string& s1, const std::string& s2, 
                         int m, int n, int startI, int startJ, 
                         const DPMatrix& scores);

    void printDPMatrix(const std::string& s1, const std::string& s2, int m, int n, const DPMatrix& score);

    void printDPAlign(const std::string& s1, const std::string& s2, 
                      int m, int n, int startI, int startJ, 
                      const DPMatrix& scores);    

};

#endif
