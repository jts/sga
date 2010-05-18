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
#include "CompleteOverlapSet.h"

// 
CompleteOverlapSet::CompleteOverlapSet(const Vertex* pVertex, double maxER, int minLength) : m_pX(pVertex), m_maxER(maxER), m_minLength(minLength)
{
    m_cost = 0;
    iterativeConstruct();
    //constructMap();
}

// Perform a breadth-first search of the graph, accumulating all the valid
// overlaps of reads to m_pX
void CompleteOverlapSet::iterativeConstruct()
{
    m_cost = 0;
    ExplorePriorityQueue queue;

    // We store the overlaps that do not meet the overlap parameters
    // in this map. The neighbors of these vertices might have valid
    // overlaps to m_pX so they must still be explored.
    SGAlgorithms::EdgeDescOverlapMap exclusionMap;

    // Add first-order overlaps
    EdgePtrVec edges = m_pX->getEdges();
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        EdgeDesc ed = pEdge->getDesc();
        Overlap ovr = pEdge->getOverlap();

        double error_rate = SGAlgorithms::calcErrorRate(m_pX, pEdge->getEnd(), ovr);
        int overlapLen = ovr.getOverlapLength(0);
        if(isErrorRateAcceptable(error_rate, m_maxER) && overlapLen >= m_minLength)
        {
            m_overlapMap.insert(std::make_pair(ed, ovr));
            queue.push(ExploreElement(ed, ovr));
        }
        else
        {
            exclusionMap.insert(std::make_pair(ed, ovr));
            queue.push(ExploreElement(ed, ovr));
        }
    }

    while(!queue.empty())
    {
        ExploreElement ee = queue.top();
        queue.pop();

        EdgeDesc& edXY = ee.ed;
        Overlap& ovrXY = ee.ovr;
        
        // Add the overlaps of Y
        Vertex* pY = edXY.pVertex;
        EdgePtrVec neighborEdges = pY->getEdges();

        for(size_t i = 0; i < neighborEdges.size(); ++i)
        {
            Edge* pEdgeYZ = neighborEdges[i];
            Vertex* pZ = pEdgeYZ->getEnd();
            if(pZ == m_pX)
                continue;

            Overlap ovrYZ = pEdgeYZ->getOverlap();
            // Check that this vertex actually overlaps pX
            if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
            {
                ++m_cost;
                Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
                
                EdgeDesc edXZ = SGAlgorithms::overlapToEdgeDesc(pZ, ovrXZ);

                double error_rate = SGAlgorithms::calcErrorRate(m_pX, pZ, ovrXZ);
                int overlapLen = ovrXZ.getOverlapLength(0);

                if(isErrorRateAcceptable(error_rate, m_maxER) && overlapLen >= m_minLength)
                {
                    // Overlap is valid
                    SGAlgorithms::EdgeDescOverlapMap::iterator findIter = m_overlapMap.find(edXZ);
                    
                    if(findIter == m_overlapMap.end())
                    {
                        m_overlapMap.insert(std::make_pair(edXZ, ovrXZ));
                        queue.push(ExploreElement(edXZ, ovrXZ));
                    }
                    else if(ovrXZ.getOverlapLength(0) > findIter->second.getOverlapLength(0))
                    {
                        findIter->second = ovrXZ;
                        queue.push(ExploreElement(edXZ, ovrXZ));
                    }   
                }
                else
                {
                    // Overlap is invalid
                    SGAlgorithms::EdgeDescOverlapMap::iterator findIter = exclusionMap.find(edXZ);
                    if(findIter == exclusionMap.end())
                    {
                        exclusionMap.insert(std::make_pair(edXZ, ovrXZ));
                        queue.push(ExploreElement(edXZ, ovrXZ));
                    }
                    else if(overlapLen > findIter->second.getOverlapLength(0))
                    {
                        findIter->second = ovrXZ;
                        queue.push(ExploreElement(edXZ, ovrXZ));
                    }
                }
            }
        }
    }
}

//
void CompleteOverlapSet::constructMap()
{
    EdgePtrVec edges = m_pX->getEdges();

    // Add the primary overlaps to the map, and all the nodes reachable from the primaries
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        EdgeDesc ed = pEdge->getDesc();
        Overlap ovr = pEdge->getOverlap();

        if(ovr.isContainment())
            continue;
        m_overlapMap.insert(std::make_pair(ed, ovr));

        // Recursively add neighbors
        recursiveConstruct(ed, ovr, 1, ovr.getOverlapLength(0));
    }
}

// Recursively add overlaps to pX inferred from the edges of pY to outMap
void CompleteOverlapSet::recursiveConstruct(const EdgeDesc& edXY, const Overlap& ovrXY, int depth, int distance)
{
    std::cout << "depth: " << depth << " " << distance << " " << ovrXY << "\n";
    Vertex* pY = edXY.pVertex;

    EdgePtrVec neighborEdges = pY->getEdges();

    for(size_t i = 0; i < neighborEdges.size(); ++i)
    {
        Edge* pEdgeYZ = neighborEdges[i];
        Vertex* pZ = pEdgeYZ->getEnd();
        if(pZ != m_pX)
        {
            Overlap ovrYZ = pEdgeYZ->getOverlap();
            if(ovrYZ.isContainment())
                continue;
            // Check that this vertex actually overlaps pX
            if(SGAlgorithms::hasTransitiveOverlap(ovrXY, ovrYZ))
            {
                Overlap ovrXZ = SGAlgorithms::inferTransitiveOverlap(ovrXY, ovrYZ);
                EdgeDesc edXZ = SGAlgorithms::overlapToEdgeDesc(pZ, ovrXZ);
                if(edXZ.dir != edXY.dir || ovrXZ.isContainment())
                    continue;

                double error_rate = SGAlgorithms::calcErrorRate(m_pX, pZ, ovrXZ);
                int overlapLen = ovrXZ.getOverlapLength(0);

                if(overlapLen > distance)
                {
                    /*
                    std::cout << "BACK OVERLAP! d: " << distance << " ol: " << overlapLen << "\n";
                    std::cout << "ovrXY " << ovrXY << "\n";
                    std::cout << "ovrYZ " << ovrYZ << "\n";
                    std::cout << "ovrXZ " << ovrXZ << "\n";
                    */
                    continue;
                }
                
                std::cout << "ovrXY " << ovrXY << "\n";
                std::cout << "ovrYZ " << ovrYZ << "\n";
                std::cout << "ovrXZ " << ovrXZ << "\n";
                if(isErrorRateAcceptable(error_rate, m_maxER) && overlapLen >= m_minLength)
                {
                    SGAlgorithms::EdgeDescOverlapMap::iterator findIter = m_overlapMap.find(edXZ);
                    
                    if(findIter == m_overlapMap.end())
                    {
                        m_overlapMap.insert(std::make_pair(edXZ, ovrXZ));
                        recursiveConstruct(edXZ, ovrXZ, depth + 1, ovrXZ.getOverlapLength(0));
                    }
                    else if(ovrXZ.getOverlapLength(0) > findIter->second.getOverlapLength(0))
                    {
                        findIter->second = ovrXZ;
                        recursiveConstruct(edXZ, ovrXZ, depth + 1, ovrXZ.getOverlapLength(0));
                    }   
                }
            }
        }
    }
}

// Compare the actual edges of m_pX with the edges in the overlap map
void CompleteOverlapSet::getDiffMap(SGAlgorithms::EdgeDescOverlapMap& missingMap, SGAlgorithms::EdgeDescOverlapMap& extraMap)
{
    missingMap = m_overlapMap;
    EdgePtrVec edges = m_pX->getEdges();

    for(size_t i = 0; i < edges.size(); ++i)
    {
        EdgeDesc ed = edges[i]->getDesc();
        SGAlgorithms::EdgeDescOverlapMap::iterator iter = missingMap.find(edges[i]->getDesc());
        if(iter != missingMap.end())
        {
            missingMap.erase(iter);
        }
        else
        {
            Overlap ovr = edges[i]->getOverlap();
            extraMap.insert(std::make_pair(ed, ovr));
        }
    }
}

// Reset the overlap parameters and filter edges
void CompleteOverlapSet::resetParameters(double maxER, int minLength)
{
    m_maxER = maxER;
    m_minLength = minLength;

    // Filter out overlaps no longer meeting the criteria
    SGAlgorithms::EdgeDescOverlapMap::iterator iter = m_overlapMap.begin();
    while(iter != m_overlapMap.end())
    {
        double er = SGAlgorithms::calcErrorRate(m_pX, iter->first.pVertex, iter->second);
        bool acceptER = isErrorRateAcceptable(er, m_maxER);
        bool acceptLen = iter->second.getOverlapLength(0) >= m_minLength;
        
        if(!acceptER || !acceptLen)
        {
            //std::cout << "Filtering " << iter->second << " flags: " << acceptER << acceptLen << "\n";
            //std::cout << "er: " << er << " len: " << iter->second.getOverlapLength(0) << "\n";
            m_overlapMap.erase(iter++);
        }
        else
            ++iter;
    }
}

// Remove all overlaps to a particular vertex
void CompleteOverlapSet::removeOverlapsTo(Vertex* pRemove)
{
    EdgeDesc ed;
    ed.pVertex = pRemove;

    for(int i = 0; i != ED_COUNT; ++i)
    {
        ed.dir = (EdgeDir)i;
        ed.comp = EC_SAME;
        m_overlapMap.erase(ed);
        ed.comp = EC_REVERSE;
        m_overlapMap.erase(ed);
    }
}

// Remove all the transitive and containment relationships, leaving only the irreducible
void CompleteOverlapSet::computeIrreducible(SGAlgorithms::EdgeDescOverlapMap* pTransitive, SGAlgorithms::EdgeDescOverlapMap* pContainments)
{
    SGAlgorithms::EdgeDescOverlapMap temp;
    partitionOverlaps(&temp, pTransitive, pContainments);
    m_overlapMap = temp;
}

// Partition the OverlapMap into edges that are containments, irreducible and transitive
// If the pointer for an output map is NULL, simply discard the edges
void CompleteOverlapSet::partitionOverlaps(SGAlgorithms::EdgeDescOverlapMap* pIrreducible, 
                                           SGAlgorithms::EdgeDescOverlapMap* pTransitive, 
                                           SGAlgorithms::EdgeDescOverlapMap* pContainment) const
{
    if(pIrreducible == NULL && pTransitive == NULL && pContainment == NULL)
        return; // Nothing to do
    m_cost = 0;
    SGAlgorithms::EdgeDescOverlapMap workingMap;

    // Stage 1, remove containments
    for(SGAlgorithms::EdgeDescOverlapMap::const_iterator iter = m_overlapMap.begin();
                                                   iter != m_overlapMap.end(); ++iter)
    {
        if(iter->second.isContainment())
        {
            if(pContainment != NULL)
                pContainment->insert(*iter);
        }
        else
        {
            workingMap.insert(*iter);
        }
    }

    // Stage 2, remove transitive edges
    SGAlgorithms::EDOPairQueue overlapQueue;
    for(SGAlgorithms::EdgeDescOverlapMap::iterator iter = workingMap.begin();
        iter != workingMap.end(); ++iter)
    {
        overlapQueue.push(*iter);
    }

    // Traverse the list of overlaps in order of length and move elements from
    // the irreducible map to the transitive map
    while(!overlapQueue.empty())
    {
        SGAlgorithms::EdgeDescOverlapPair edoPair = overlapQueue.top();
        overlapQueue.pop();

        EdgeDesc& edXY = edoPair.first;
        Overlap& ovrXY = edoPair.second;

        assert(!ovrXY.match.isContainment());

        SGAlgorithms::EdgeDescOverlapMap::iterator iter = workingMap.begin();
        while(iter != workingMap.end())
        {
            bool move = false;
            const EdgeDesc& edXZ = iter->first;
            const Overlap& ovrXZ = iter->second;
            m_cost++;
            // Four conditions must be met to mark an edge X->Z transitive through X->Y
            // 1) The overlaps must be in the same direction
            // 2) The overlap X->Y must be strictly longer than X->Z
            // 3) The overlap between Y->Z must not be a containment
            // 4) The overlap between Y->Z must be within the error and length thresholds
            if(!(edXZ == edXY) && edXY.dir == edXZ.dir && ovrXY.getOverlapLength(0) > ovrXZ.getOverlapLength(0))
            {
                // Infer the YZ overlap
                Overlap ovrYX = ovrXY;
                ovrYX.swap();

                Overlap ovrYZ = SGAlgorithms::inferTransitiveOverlap(ovrYX, ovrXZ);
                if(!ovrYZ.match.isContainment())
                {  
                    // Compute the error rate between the sequences
                    double error_rate = SGAlgorithms::calcErrorRate(edXY.pVertex, edXZ.pVertex, ovrYZ);
                    
                    //std::cout << "\tOVRXY: " << ovrXY << "\n";
                    //std::cout << "\tOVRXZ: " << ovrXZ << "\n";
                    //std::cout << "\tOVRYZ: " << ovrYZ << " er: " << error_rate << "\n";
                    
                    if(isErrorRateAcceptable(error_rate, m_maxER) && ovrYZ.getOverlapLength(0) >= m_minLength)
                        move = true;
                }
            }
            
            if(move)
            {
                //std::cout << "Marking overlap: " << iter->second << " as trans via " << ovrXY << "\n";
                if(pTransitive != NULL)
                    pTransitive->insert(*iter);
                workingMap.erase(iter++);
            }
            else
            {
                ++iter;
            }
        }
    }

    *pIrreducible = workingMap;
}
