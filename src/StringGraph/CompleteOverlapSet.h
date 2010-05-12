//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// CompleteOverlapSet - A collection of all
// valid overlaps for a given vertex. This is
// inferred from the structure of a transitive reduced
// graph
//
#ifndef COMPLETEOVERLAPSET_H
#define COMPLETEOVERLAPSET_H

#include "Bigraph.h"
#include "SGAlgorithms.h"

class CompleteOverlapSet
{
    public:

        CompleteOverlapSet(Vertex* pVertex, double maxER, int minLength);

        void getDiffMap(SGAlgorithms::EdgeDescOverlapMap& missingMap, SGAlgorithms::EdgeDescOverlapMap& extraMap);
        void removeOverlapsTo(Vertex* pRemove);
        void removeTransitiveOverlaps();
        void resetParameters(double maxER, int minLength);

    private:

        // functions
        void constructMap();
        void recursiveConstruct(const EdgeDesc& edXY, const Overlap& ovrXY);

        // data
        SGAlgorithms::EdgeDescOverlapMap m_overlapMap;

        // the vertex this set is centered on
        Vertex* m_pX;

        double m_maxER;
        int m_minLength;

};

#endif
