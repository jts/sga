//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGSearch - Algorithms and data structures
// for searching a string graph
//
#include "SGSearch.h"
#include <queue>

//
SGWalk::SGWalk(const Vertex* pStartVertex) : m_pStartVertex(pStartVertex), m_extensionDistance(0)
{

}

//
void SGWalk::addEdge(Edge* pEdge)
{
    m_edges.push_back(pEdge);
    m_extensionDistance += pEdge->getSeqLen();
}

// Construct the extension string corresponding to the path
std::string SGWalk::getString(SGWalkType type) const
{
    std::string out;

    if(type == SGWT_START_TO_END)
    {
        out.append(m_pStartVertex->getSeq().toString());
    }


    // Determine if the string should go to the end of the last vertex
    // in the path
    size_t stop = m_edges.size();

    // The first edge is always in correct frame of reference 
    // so the comp is EC_SAME. This variable tracks where the 
    // string that is being added is different from the starting sequence
    // and needs to be flipped
    EdgeComp currComp = EC_SAME;

    // If the walk direction is antisense, we reverse every component and then
    // reverse the entire string to generate the final string
    bool reverseAll = !m_edges.empty() && m_edges[0]->getDir() == ED_ANTISENSE;
    if(reverseAll)
        out = reverse(out);

    for(size_t i = 0; i < stop; ++i)
    {
        Edge* pYZ = m_edges[i];
        
        // Append in the extension string
        std::string edge_str = pYZ->getLabel();
        assert(edge_str.size() != 0);
        if(currComp == EC_REVERSE)
            edge_str = reverseComplement(edge_str);

        if(reverseAll)
            edge_str = reverse(edge_str);

        // Calculate the next comp, between X and Z
        EdgeComp ecYZ = pYZ->getComp();
        EdgeComp ecXZ;
        if(ecYZ == EC_SAME)
            ecXZ = currComp;
        else
            ecXZ = !currComp;

        out.append(edge_str);
        currComp = ecXZ;
    }

    if(reverseAll)
        out = reverse(out);
    return out;
}

//
int SGWalk::getExtensionDistance() const
{
    return m_extensionDistance;
}

// This is equivalent to the extension distance
int SGWalk::getEndToEndDistance() const
{
    return m_extensionDistance;
}

// 
int SGWalk::getStartToEndDistance() const
{
    return m_pStartVertex->getSeqLen() + getEndToEndDistance();
}

// Returns the distance from the end of pStart to the beginning of the last vertex in the path
// This can be negative if they overlap
int SGWalk::getEndToStartDistance() const
{
    if(m_edges.empty())
        return 0;

    int len_x = m_pStartVertex->getSeqLen();
    int len_y = getLastEdge()->getEnd()->getSeqLen();
    return getStartToEndDistance() - (len_x + len_y);
}

//
Edge* SGWalk::getLastEdge() const
{
    assert(!m_edges.empty());
    return m_edges.back();
}

// 
void SGWalk::print() const
{
    std::cout << "Walk start: " << m_pStartVertex->getID() << "\nWalk: ";
    const Vertex* pLast = m_pStartVertex;
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
    {
        //std::cout << *(*iter) << " ";
        std::cout << (*iter)->getStartID() << " -- " << (*iter)->getEndID() << "," << (*iter)->getDir() << "," << (*iter)->getComp() << "\t";
        assert((*iter)->getStart() == pLast);
        pLast = (*iter)->getEnd();
    }
    std::cout << "\n";
}   

// Find all the walks between pX and pY that are within maxDistance
void SGSearch::findWalks(const Vertex* pX, const Vertex* pY, EdgeDir initialDir,
                         int maxDistance, size_t maxQueue, SGWalkVector& outWalks)
{
    // Create the initial path nodes
    typedef std::queue<SGWalk> WalkQueue;
    EdgeDir dir = initialDir;
    WalkQueue queue;
    EdgePtrVec edges = pX->getEdges(dir);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        assert(!pEdge->getOverlap().isContainment());

        SGWalk walk(pX);
        walk.addEdge(pEdge);
        queue.push(walk);
    }

    while(!queue.empty())
    {
        if(queue.size() > maxQueue)
        {
            // Give up the search if there are too many possible paths to continue
            std::cerr << "Warning queue too large " << queue.size() << "\n";
            outWalks.clear();
            return;
        }

        bool bPop = false;
        
        // Process the current walk
        // This occurs in a lower scope so we can safely use a reference
        // to the top element before popping it off without having the 
        // reference hang around.
        {
            SGWalk& currWalk = queue.front();
            Edge* pWZ = currWalk.getLastEdge(); 
            Vertex* pZ = pWZ->getEnd();
           
            // Check if we have found pY or exceeded the distance
            if(pZ == pY)
            {
                outWalks.push_back(currWalk);
                bPop = true;
            }
            else if(currWalk.getExtensionDistance() > maxDistance)
            {
                bPop = true;
            }
            else
            {
                // Path is still valid, continue
                EdgeDir continueDir = pWZ->getTransitiveDir();
                EdgePtrVec zEdges = pZ->getEdges(continueDir);

                if(!zEdges.empty())
                {
                    // If there is one extension, just add it to the current walk
                    if(zEdges.size() == 1)
                    {
                        Edge* pNextEdge = zEdges[0];
                        currWalk.addEdge(pNextEdge);
                        bPop = false;
                    }
                    else
                    {
                        // Multiple extensions, remove the current walk and branch
                        bPop = true;
                        for(size_t i = 0; i < zEdges.size(); ++i)
                        {
                            SGWalk branch = currWalk;
                            Edge* pBranchEdge = zEdges[i];
                            branch.addEdge(pBranchEdge);
                            queue.push(branch);
                        }
                    }
                }
                else
                {
                    bPop = true;
                }
            }
        }

        if(bPop)
            queue.pop();
    }
}

