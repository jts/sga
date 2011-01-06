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
SGWalk::SGWalk(Vertex* pStartVertex, bool bIndexWalk) : m_pStartVertex(pStartVertex), m_extensionDistance(0), m_extensionFinished(false)
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
    m_extensionFinished = other.m_extensionFinished;
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
void SGWalk::setFinished(bool b)
{
    m_extensionFinished = b;
}

bool SGWalk::isFinished() const
{
    return m_extensionFinished;
}

//
void SGWalk::popLast()
{
    m_edges.pop_back();
}

// 
Vertex* SGWalk::getStartVertex() const
{
    return m_pStartVertex;
}

//
size_t SGWalk::getNumVertices() const
{
    return getNumEdges() + 1;
}

//
size_t SGWalk::getNumEdges() const
{
    return m_edges.size();
}

bool SGWalk::isIndexed() const
{
    return m_pWalkIndex != NULL;
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
VertexPtrVec SGWalk::getVertices() const
{
    VertexPtrVec out;
    out.push_back(m_pStartVertex);
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        out.push_back((*iter)->getEnd());
    return out;
}

// Return a string describing this path
std::string SGWalk::pathSignature() const
{
    std::stringstream ss;
    ss << m_pStartVertex->getID() << ",";
    for(EdgePtrVec::const_iterator iter = m_edges.begin(); iter != m_edges.end(); ++iter)
        ss << (*iter)->getEndID() << ", ";
    return ss.str();
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

// 
void SGWalk::printSimple() const
{
    std::cout << pathSignature() << "\n";
}

// Find all the walks between pX and pY that are within maxDistance
void SGSearch::findWalks(Vertex* pX, Vertex* pY, EdgeDir initialDir,
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

// Search the graph for a set of walks that represent alternate
// versions of the same sequence. Theese walks are found by searching
// the graph for a set of walks that start/end at a common vertex and cover
// every vertex inbetween the start/end points. Additionally, the internal
// vertices cannot have an edge to any vertex that is not in the set of vertices
// indicated by the walks.
// If these two conditions are met, we can retain one of the walks and cleanly remove
// the others.
void SGSearch::findVariantWalks(Vertex* pX, 
                                EdgeDir initialDir, 
                                int maxDistance,
                                size_t maxWalks, 
                                SGWalkVector& outWalks)
{
    findCollapsedWalks(pX, initialDir, maxDistance, maxWalks, outWalks);

    if(outWalks.empty())
        return;
    assert(outWalks.size() > 1);

    // Validate that any of the returned walks can be removed cleanly.
    // This means that for all the internal vertices on each walk (between
    // the common endpoints) the only links are to other vertices in the set

    // Construct the set of vertex IDs
    std::set<Vertex*> completeVertexSet;
    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        VertexPtrVec verts = outWalks[i].getVertices();
        for(size_t j = 0; j < verts.size(); ++j)
        {
            completeVertexSet.insert(verts[j]);
        }
    }

    // Check that each vertex in the internal nodes only has
    // edges to the vertices in the set

    bool cleanlyRemovable = true;

    // Ensure that all the vertices linked to the start vertex
    // in the specified dir are present in the set.
    EdgePtrVec epv = pX->getEdges(initialDir);
    cleanlyRemovable = checkEndpointsInSet(epv, completeVertexSet);

    // Ensure that all the vertex linked to the last vertex
    // in the incoming direction are preset
    Edge* pLastEdge = outWalks.front().getLastEdge();
    Vertex* pLastVertex = pLastEdge->getEnd();
    EdgeDir lastDir = pLastEdge->getTwinDir();
    
    epv = pLastVertex->getEdges(lastDir);

    cleanlyRemovable = cleanlyRemovable && checkEndpointsInSet(epv, completeVertexSet);

    // Check that each vertex connected to an interval vertex is also present
    for(std::set<Vertex*>::iterator iter = completeVertexSet.begin(); iter != completeVertexSet.end(); ++iter)
    {
        Vertex* pY = *iter;
        if(pY == pX || pY == pLastVertex)
            continue;
        epv = pY->getEdges();
        cleanlyRemovable = cleanlyRemovable && checkEndpointsInSet(epv, completeVertexSet);
    }

    if(!cleanlyRemovable)
    {
        outWalks.clear();
    }
    return;
}

// Check that all the endpoints of the edges in the edge pointer vector
// are members of the set
bool SGSearch::checkEndpointsInSet(EdgePtrVec& epv, std::set<Vertex*>& vertexSet)
{
    for(size_t i = 0; i < epv.size(); ++i)
    {
        if(vertexSet.find(epv[i]->getEnd()) == vertexSet.end())
            return false;
    }
    return true;
}

// Return a set of walks that all start from pX and join together at some later vertex
// If no such walk exists, an empty set is returned
void SGSearch::findCollapsedWalks(Vertex* pX, EdgeDir initialDir, 
                                  int /*maxDistance*/, size_t maxQueue, 
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

        bool allFinished = true;
        for(size_t i = 0; i < queue.size(); ++i)
        {
            // Check if this path can no longer be extended
            allFinished = allFinished && queue[i].isFinished();

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

                // After truncation, redundant paths may exist if the common node had a branch.
                // We remove the redundant paths by inserting the paths into a map based on their signature string
                std::map<std::string, SGWalk*> nonRedundant;
                for(size_t i = 0; i < queue.size(); ++i)
                {
                    // Truncate the path at the common vertex
                    queue[i].truncate(iLastID);
                    SGWalk* pWalk = &queue[i];
                    nonRedundant.insert(std::make_pair(queue[i].pathSignature(), pWalk));
                }

                // Write out the paths
                for(std::map<std::string, SGWalk*>::iterator iter = nonRedundant.begin();
                    iter != nonRedundant.end(); ++iter)
                {
                    outWalks.push_back(*iter->second);
                }
                return;
            }

            if(allFinished)
            {
                // no more extension can occur, return an empty set
                outWalks.clear();
                return;
            }
        }

        if(queue.size() > maxQueue)
        {
            // Give up the search if there are too many possible paths to continue
            outWalks.clear();
            return;
        }

        // Extend each walk in the queue
        // If the walk branches, the new branch is added to the incoming queue
        // which is merged into the full queue after the extension round
        WalkQueue incoming;
        for(size_t i = 0; i < queue.size(); ++i)
        {
            SGWalk& currWalk = queue[i];
            Edge* pWZ = currWalk.getLastEdge(); 
            Vertex* pZ = pWZ->getEnd();
            extendWalk(pZ, pWZ->getTransitiveDir(), currWalk, incoming);
        }

        // Copy any new branches into the queue
        queue.insert(queue.end(), incoming.begin(), incoming.end());
    }

    (void)outWalks;
}

// Count the number of reads that span the junction described by edge XY
// X --------------
// Y        ------------
// Z            -----------
// W                ----------
// In this case Z spans the junction but W does not. 
int SGSearch::countSpanningCoverage(Edge* pXY, size_t maxQueue)
{
    Vertex* pX = pXY->getStart();

    SGWalk walk(pX, false);
    walk.addEdge(pXY);
    
    WalkQueue queue;
    queue.push_back(walk);

    // Create the initial queue
    SGWalkVector outWalks;

    //
    while(queue.size() > 0)
    {
        if(queue.size() > maxQueue)
        {
            // Give up the search if there are too many possible paths to continue
            return -1;
        }

        bool bPop = false;
        {
            SGWalk& currWalk = queue.front();
            Edge* pWZ = currWalk.getLastEdge(); 
            Vertex* pZ = pWZ->getEnd();

            // Calculate the length of the overlap of pZ on pX. If it is <= 0
            // pZ does not overlap pX and will be removed from the walk and the walk
            // is not processed further
            int extensionDistance = currWalk.getExtensionDistance();
            int overlap = pZ->getSeqLen() - extensionDistance;
            if(overlap <= 0)
            {
                //std::cout << "Too far: " << currWalk.getExtensionDistance() << "\n";
                bPop = true;
                currWalk.popLast();
                outWalks.push_back(currWalk);
            }
            else
            {
                // Continue walk
                bPop = !extendWalk(pZ, pWZ->getTransitiveDir(), currWalk, queue);
            }
        }

        if(bPop)
            queue.pop_front();
    }

    // The outwalks may have redundant sequences
    // Make a set to calculate the coverage by unique vertices
    std::set<VertexID> vertexSet;

    for(size_t i = 0; i < outWalks.size(); ++i)
    {
        SGWalk& walk = outWalks[i];
        for(size_t j = 0; j < walk.getNumEdges(); ++j)
        {
            vertexSet.insert(walk.getEdge(j)->getEndID());
        }
    }

    return vertexSet.size();
}

//
void SGSearch::initializeWalkQueue(Vertex* pX, EdgeDir initialDir, bool bIndexWalks, WalkQueue& queue)
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
    {
        currWalk.setFinished(true);
        return false;
    }

    // If there are multiple extensions, create new branches and add them to the queue
    for(size_t i = 1; i < edges.size(); ++i)
    {
        SGWalk branch = currWalk; 
        Edge* pBranchEdge = edges[i];
        branch.addEdge(pBranchEdge);
        queue.push_back(branch);
    }

    // Extend the current walk with the first edge
    Edge* pFirstEdge = edges[0];
    currWalk.addEdge(pFirstEdge);
    return true;
}
