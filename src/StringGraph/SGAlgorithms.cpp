//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGAlgorithms - Collection of algorithms for operating on string graphs
//
#include "SGAlgorithms.h"
#include "SGUtil.h"
#include "ErrorCorrect.h"
#include "CompleteOverlapSet.h"
#include <iterator>

// Find new edges for pVertex that are required if pDeleteEdge is removed from the graph
void SGAlgorithms::remodelVertexForExcision(StringGraph* pGraph, Vertex* pVertex, Edge* pDeleteEdge)
{
    assert(pVertex == pDeleteEdge->getStart());
    double maxER = pGraph->getErrorRate();
    int minLength = pGraph->getMinOverlap();

    CompleteOverlapSet vertexOverlapSet(pVertex, maxER, minLength);
    vertexOverlapSet.removeOverlapsTo(pDeleteEdge->getEnd());
    vertexOverlapSet.removeTransitiveOverlaps();
    EdgeDescOverlapMap addMap;
    EdgeDescOverlapMap removeMap;
    vertexOverlapSet.getDiffMap(addMap, removeMap);
       
    //assert(removeMap.size() == 0);
    for(EdgeDescOverlapMap::iterator iter = addMap.begin();
        iter != addMap.end(); ++iter)
    {
        //std::cout << "Adding edge " << iter->second << "\n";
        SGUtil::createEdges(pGraph, iter->second, true);
    }
}

// Calculate the error rate between the two vertices
double SGAlgorithms::calcErrorRate(const Vertex* pX, const Vertex* pY, const Overlap& ovrXY)
{
    int num_diffs = ovrXY.match.countDifferences(pX->getSeq(), pY->getSeq());
    return static_cast<double>(num_diffs) / static_cast<double>(ovrXY.match.getMinOverlapLength());
}

// Infer an overlap from two edges
// The input edges are between X->Y Y->Z
// and the returned overlap is X->Z
Overlap SGAlgorithms::inferTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ)
{
    // Construct the match
    Match match_yx = ovrXY.match;
    match_yx.swap(); 
    Match match_yz = ovrYZ.match;

    // Infer the match_ij based match_i and match_j
    Match match_xz = Match::infer(match_yx, match_yz);
    match_xz.expand();

    // Convert the match to an overlap
    Overlap ovr(ovrXY.id[0], ovrYZ.id[1], match_xz);
    return ovr;
}

// Return a descriptor of the edge describing ovrXY
EdgeDesc SGAlgorithms::overlapToEdgeDesc(Vertex* pY, const Overlap& ovrXY)
{
    EdgeDesc edXY;
    edXY.pVertex = pY;
    edXY.comp = (ovrXY.match.isRC()) ? EC_REVERSE : EC_SAME;
    edXY.dir = ovrXY.match.coord[0].isLeftExtreme() ? ED_ANTISENSE : ED_SENSE; // X -> Y
    return edXY;
}

// Return true if XZ has an overlap
bool SGAlgorithms::hasTransitiveOverlap(const Overlap& ovrXY, const Overlap& ovrYZ)
{
    Match match_yx = ovrXY.match;
    match_yx.swap(); 
    Match match_yz = ovrYZ.match;
    return Match::doMatchesIntersect(match_yx, match_yz);
}

// Construct an extended multioverlap for a vertex
MultiOverlap SGAlgorithms::makeExtendedMultiOverlap(const StringGraph* pGraph, const Vertex* pVertex)
{
    CompleteOverlapSet overlapSet(pVertex, pGraph->getErrorRate(), 1);
    EdgeDescOverlapMap overlapMap = overlapSet.getOverlapMap();

    MultiOverlap mo(pVertex->getID(), pVertex->getSeq());
    for(EdgeDescOverlapMap::const_iterator iter = overlapMap.begin();
        iter != overlapMap.end(); ++iter)
    {
        mo.add(iter->first.pVertex->getSeq(), iter->second);
    }
    return mo;
}

//
void SGAlgorithms::makeExtendedSeqTries(const StringGraph* pGraph, const Vertex* pVertex, double p_error, SeqTrie* pLeftTrie, SeqTrie* pRightTrie)
{
    double lp = log(p_error);
    CompleteOverlapSet overlapSet(pVertex, pGraph->getErrorRate(), 1);
    EdgeDescOverlapMap overlapMap = overlapSet.getOverlapMap();

    for(EdgeDescOverlapMap::const_iterator iter = overlapMap.begin();
        iter != overlapMap.end(); ++iter)
    {
        // Coord[0] of the match is wrt pVertex, coord[1] is the other read
        std::string overlapped = iter->second.match.coord[1].getSubstring(iter->first.pVertex->getSeq());
        if(iter->second.match.isRC())
            overlapped = reverseComplement(overlapped);

        if(iter->second.match.coord[0].isRightExtreme())
        {
            overlapped = reverse(overlapped);
            pRightTrie->insert(overlapped, lp);
        }
        else if(iter->second.match.coord[0].isLeftExtreme())
        {
            pLeftTrie->insert(overlapped, lp);
        }
        else
        {
            std::cout << "FOUND CONTAINMENT: " << iter->second << "\n";
            assert(iter->second.match.coord[0].isLeftExtreme());
        }
    }        
}

//
EdgeDesc SGAlgorithms::getEdgeDescFromEdge(Edge* pEdge)
{
    return pEdge->getDesc();
}

//
EdgeDesc SGAlgorithms::getEdgeDescFromPair(const EdgeDescOverlapPair& pair)
{
    return pair.first;
}

void SGAlgorithms::printOverlapMap(const EdgeDescOverlapMap& overlapMap)
{
    for(EdgeDescOverlapMap::const_iterator iter = overlapMap.begin(); iter != overlapMap.end(); ++iter)
    {
        std::cout << "  Overlap:" << iter->second << "\n";
    }
}

