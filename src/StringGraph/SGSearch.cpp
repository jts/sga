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
SGWalk::SGWalk(const Vertex* pStartVertex, bool bIndexWalk) : m_pStartVertex(pStartVertex), m_extensionDistance(0)
{
    if(bIndexWalk)
        m_pWalkIndex = new WalkIndex;
    else
        m_pWalkIndex = NULL;
}

//
SGWalk::~SGWalk()
{
    if(m_pWalkIndex != NULL)
    {
        delete m_pWalkIndex;
        m_pWalkIndex = NULL;
    }
}

//
SGWalk::SGWalk(const SGWalk& other)
{
    m_pStartVertex = other.m_pStartVertex;
    m_edges = other.m_edges;
    m_extensionDistance = other.m_extensionDistance;
    m_pWalkIndex = NULL;

    if(other.m_pWalkIndex != NULL)
        m_pWalkIndex = new WalkIndex(*other.m_pWalkIndex);
}

//
SGWalk& SGWalk::operator=(const SGWalk& other)
{
    if(&other == this)
        return *this; // self-assign
   
    m_pStartVertex = other.m_pStartVertex;
    m_edges = other.m_edges;
    m_extensionDistance = other.m_extensionDistance;
    
    if(m_pWalkIndex != NULL)
        delete m_pWalkIndex;
    
    if(other.m_pWalkIndex != NULL)
        m_pWalkIndex = new WalkIndex(*other.m_pWalkIndex);
        
    return *this;
}

//
void SGWalk::addEdge(Edge* pEdge)
{
    m_edges.push_back(pEdge);
    m_extensionDistance += pEdge->getSeqLen();

    // Add the vertex ID to the index if necessary
    if(m_pWalkIndex != NULL)
        m_pWalkIndex->insert(m_edges.back()->getEndID());
}

//
size_t SGWalk::getNumEdges() const
{
    return m_edges.size();
}

//
bool SGWalk::containsVertex(const VertexID& id) const
{
    if(m_pWalkIndex == NULL)
        assert(false);
    return m_pWalkIndex->count(id) > 0;
}

//
void SGWalk::truncate(const VertexID& id)
{
    EdgePtrVec::iterator iter = m_edges.begin();
    while(iter != m_edges.end())
    {
        if((*iter)->getEndID() == id)
            break;
        ++iter;
    }
    assert(iter != m_edges.end());
    m_edges.erase(iter + 1, m_edges.end());
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
Edge* SGWalk::getEdge(size_t idx) const
{
    return m_edges[idx];
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
    WalkQueue queue;
    initializeWalkQueue(pX, initialDir, false, queue);

    while(!queue.empty())
    {
        if(queue.size() > maxQueue)
        {
            // Give up the search if there are too many possible paths to continue
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
                bPop = !extendWalk(pZ, pWZ->getTransitiveDir(), currWalk, queue);
            }
        }

        if(bPop)
            queue.pop_front();
    }
}

// Return a set of walks that all start from pX and join together at some later vertex
// If no such walk exists, an empty set is returned
void SGSearch::findCollapsedWalks(const Vertex* pX, EdgeDir initialDir, 
                                  int maxDistance, size_t maxQueue, 
                                  SGWalkVector& outWalks)
{
    // Create the initial path nodes
    WalkQueue queue;
    initializeWalkQueue(pX, initialDir, true, queue);

    //
    while(queue.size() > 1)
    {
        // Check the last element of each walk in the queue to see if all
        // the walks share a common vertex
        for(size_t i = 0; i < queue.size(); ++i)
        {
            VertexID iLastID = queue[i].getLastEdge()->getEndID();
            bool isCommonVertex = true;
            for(size_t j = 0; j < queue.size(); ++j)
            {
                if(j == i)
                    continue;
                if(!queue[j].containsVertex(iLastID))
                {
                    isCommonVertex = false;
                    break;
                }
            }

            if(isCommonVertex)
            {
                // This vertex is common between all walks, return the found walks in outWalks
                //std::cout << "Vertex " << iLastID << " is common to all walks\n";
                for(size_t i = 0; i < queue.size(); ++i)
                {
                    // Truncate the path at the common vertex
                    queue[i].truncate(iLastID);
                    outWalks.push_back(queue[i]);
                }
                return;
            }
        }

        if(queue.size() > maxQueue)
        {
            // Give up the search if there are too many possible paths to continue
            outWalks.clear();
            return;
        }

        bool bPop = false;
        {
            SGWalk& currWalk = queue.front();
            Edge* pWZ = currWalk.getLastEdge(); 
            Vertex* pZ = pWZ->getEnd();
            if(currWalk.getExtensionDistance() > maxDistance)
            {
                //std::cout << "Too far: " << currWalk.getExtensionDistance() << "\n";
                bPop = true;
            }
            else
            {
                bPop = !extendWalk(pZ, pWZ->getTransitiveDir(), currWalk, queue);
            }
        }

        if(bPop)
            queue.pop_front();
    }

    (void)outWalks;
}

//
void SGSearch::initializeWalkQueue(const Vertex* pX, EdgeDir initialDir, bool bIndexWalks, WalkQueue& queue)
{
    EdgePtrVec edges = pX->getEdges(initialDir);
    for(size_t i = 0; i < edges.size(); ++i)
    {
        Edge* pEdge = edges[i];
        assert(!pEdge->getOverlap().isContainment());

        SGWalk walk(pX, bIndexWalks);
        walk.addEdge(pEdge);
        queue.push_back(walk);
    }
}

// Extend the walk by addding the neighbors of the last vertex (pX) to the walk
// Returns true if the branch was extended
bool SGSearch::extendWalk(const Vertex* pX, EdgeDir dir, SGWalk& currWalk, WalkQueue& queue)
{
    EdgePtrVec edges = pX->getEdges(dir);

    if(edges.empty())
        return false;
    
    // If there are multiple extensions, create new branches and add them to the queue
    for(size_t i = 1; i < edges.size(); ++i)
    {
        //std::cout << "Branching from " << pX->getID() << " to " << edges[i]->getEndID() << "\n";
        SGWalk branch = currWalk; 
        Edge* pBranchEdge = edges[i];
        branch.addEdge(pBranchEdge);
        queue.push_back(branch);
    }

    // Extend the current walk with the first edge
    //std::cout << "Extending from " << pX->getID() << " to " << edges[0]->getEndID() << "\n";
    Edge* pFirstEdge = edges[0];
    currWalk.addEdge(pFirstEdge);
    return true;
}
