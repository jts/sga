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

// Structure used for iterative exploration of the graph
// The exploration starts at some vertex X, each element
// holds a possible overlap between X and another vertex Y
struct ExploreElement
{
    ExploreElement(const EdgeDesc& e, const Overlap& o) : ed(e), ovr(o) {}
    EdgeDesc ed;
    Overlap ovr;
};

// Comparison operator used to compare ExploreElements
// by the length of the overlap on vertex X
struct CompareExploreElemOverlapLength
{
    bool operator()(const ExploreElement& elemXY, const ExploreElement& elemXZ)
    {
        return elemXY.ovr.match.coord[0].length() < elemXZ.ovr.match.coord[0].length();
    }
};

//
typedef std::priority_queue<ExploreElement, 
                            std::vector<ExploreElement>, 
                            CompareExploreElemOverlapLength> ExplorePriorityQueue;


class CompleteOverlapSet
{
    public:

        CompleteOverlapSet(const Vertex* pVertex, double maxER, int minLength);

        void getDiffMap(SGAlgorithms::EdgeDescOverlapMap& missingMap, SGAlgorithms::EdgeDescOverlapMap& extraMap);
        void removeOverlapsTo(Vertex* pRemove);
        void removeTransitiveOverlaps();
        void resetParameters(double maxER, int minLength);
        SGAlgorithms::EdgeDescOverlapMap getOverlapMap() const { return m_overlapMap; }

    private:

        // functions
        void iterativeConstruct();

        // data
        SGAlgorithms::EdgeDescOverlapMap m_overlapMap;

        // the vertex this set is centered on
        const Vertex* m_pX;

        double m_maxER;
        int m_minLength;

};

#endif
