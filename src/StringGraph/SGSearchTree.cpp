//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGSearchTree - Data structure used for implementing
// graph search algorithms
//
#include "SGSearchTree.h"

//
SGSearchNode::SGSearchNode(Vertex* pVertex,
                           EdgeDir expandDir,
                           SGSearchNode* pParent,
                           Edge* pEdgeFromParent) : m_pVertex(pVertex),
                                                    m_expandDir(expandDir),
                                                    m_pParent(pParent), 
                                                    m_pEdgeFromParent(pEdgeFromParent),
                                                    m_numChildren(0)
{
    // Set the extension distance
    if(m_pParent == NULL)
    {
        // this is the root node with extension distance 0
        m_extensionDistance = 0;
    }
    else
    {
        assert(m_pEdgeFromParent != NULL);
        m_extensionDistance = m_pParent->m_extensionDistance + m_pEdgeFromParent->getSeqLen();
    }
}

// Delete this node and decrement the number of children
// in the parent node. All children of a node must
// be deleted before the parent
SGSearchNode::~SGSearchNode()
{
    assert(m_numChildren == 0);
    if(m_pParent != NULL)
        m_pParent->decrementChildren();
}

//
void SGSearchNode::decrementChildren()
{
    assert(m_numChildren != 0);
    m_numChildren -= 1;
}

// creates nodes for the children of this node
// and place pointers to them in the queue.
// Returns the number of nodes created
int SGSearchNode::createChildren(SGSearchNodePtrDeque& outDeque)
{
    assert(m_numChildren == 0);

    EdgePtrVec edges = m_pVertex->getEdges(m_expandDir);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        EdgeDir childExpandDir = !edges[i]->getTwin()->getDir();
        SGSearchNode* pNode = new SGSearchNode(edges[i]->getEnd(), childExpandDir, this, edges[i]);
        outDeque.push_back(pNode);
        m_numChildren += 1;
    }
    return edges.size();
}

//
// SGSearchTree
//
SGSearchTree::SGSearchTree(Vertex* pStartVertex, 
                           Vertex* pEndVertex, 
                           EdgeDir searchDir,
                           size_t distanceLimit,
                           size_t nodeLimit) : m_pGoalVertex(pEndVertex),
                                            m_distanceLimit(distanceLimit),
                                            m_nodeLimit(nodeLimit),
                                            m_searchAborted(false)
{
    // Create the root node of the search tree
    m_pRootNode = new SGSearchNode(pStartVertex, searchDir, NULL, NULL);

    // add the root to the expand queue
    m_expandQueue.push_back(m_pRootNode);

    m_totalNodes = 1;
}

SGSearchTree::~SGSearchTree()
{
    // Delete the tree
    // We delete each leaf and recurse up the tree iteratively deleting
    // parents with a single child node. This ensure that each parent is
    // deleted after all its children
    SGSearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);
 
    size_t totalDeleted = 0;
    for(SGSearchNodePtrDeque::iterator iter = m_expandQueue.begin(); 
                                       iter != m_expandQueue.end();
                                       ++iter)    
    {
        SGSearchNode* pCurr = *iter;

        // loop invariant: pCurr is a deletable node
        // the loop stops when the parent is NULL or has 
        // a child other than pCurr
        do
        {
            assert(pCurr->getNumChildren() == 0);
            SGSearchNode* pNext = pCurr->getParent();
            
            delete pCurr; // decrements pNext's child count
            totalDeleted += 1;

            pCurr = pNext;
        } while(pCurr && pCurr->getNumChildren() == 0);
    }

    assert(totalDeleted == m_totalNodes);
}

// Perform one step of the BFS
bool SGSearchTree::stepOnce()
{
    if(m_expandQueue.empty())
        return false;

    if(m_totalNodes > m_nodeLimit)
    {
        // Move all nodes in the expand queue to the done queue
        m_doneQueue.insert(m_doneQueue.end(), m_expandQueue.begin(), m_expandQueue.end());
        m_expandQueue.clear();

        // Set a flag indicating the search was aborted
        m_searchAborted = true;
        return false;
    }

    // Iterate over the expand queue. If the path to the node
    // is outside the depth limit, move that node to the done queue. It cannot
    // yield a valid path to the goal. Otherwise, add the children of the node 
    // to the incoming queue
    SGSearchNodePtrDeque incomingQueue;
    while(!m_expandQueue.empty())
    {
        SGSearchNode* pNode = m_expandQueue.front();
        m_expandQueue.pop_front();

        if(pNode->getVertex() == m_pGoalVertex)
        {
            // This node represents the goal, add it to the goal queue
            m_goalQueue.push_back(pNode);
            continue;
        }

        if(pNode->getDistance() > m_distanceLimit)
        {
            // Path to this node is too long, expand it no further
            m_doneQueue.push_back(pNode);
        }
        else
        {
            // Add the children of this node to the queue
            int numCreated = pNode->createChildren(incomingQueue);
            m_totalNodes += numCreated;

            if(numCreated == 0)
            {
                // No children created, add this node to the done queue
                m_doneQueue.push_back(pNode);
            }
        }
    }

    m_expandQueue = incomingQueue;
    return true;
}

// Return true if all the walks from the root converge
// to one vertex (ie if the search from pX converged
// to pY, then ALL paths from pX must go through pY).
bool SGSearchTree::hasSearchConverged(Vertex*& pConvergedVertex)
{
    // Construct a set of all the leaf nodes
    SGSearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    // Search all the tree for all the nodes in the expand queue
    for(SGSearchNodePtrDeque::iterator iter = m_expandQueue.begin(); 
                                       iter != m_expandQueue.end();
                                       ++iter)
    {
        SGSearchNode* pNode = *iter;

        bool isInAllBranches = true;
        for(SGSearchNodePtrDeque::iterator leafIter = completeLeafNodes.begin();
                                           leafIter != completeLeafNodes.end();
                                           ++leafIter)
        {
            // Search the current branch from this leaf node to the root
            bool isInBranch = searchBranchForVertex(*leafIter, pNode->getVertex());
            if(!isInBranch)
            {
                isInAllBranches = false;
                break;
            }
        }
    
        // search has converted
        if(isInAllBranches)
        {
            pConvergedVertex = pNode->getVertex();
            return true;
        }
    }

    // search has not converted
    pConvergedVertex = NULL;
    return false;
}

// Construct walks representing every path from the start node
void SGSearchTree::buildWalksToAllLeaves(SGWalkVector& outWalks)
{
    // Construct a queue with all leaf nodes in it
    SGSearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    _buildWalksToLeaves(completeLeafNodes, outWalks);
}

// Construct walks representing every path from the start vertex to the goal vertex
void SGSearchTree::buildWalksToGoal(SGWalkVector& outWalks)
{
    _buildWalksToLeaves(m_goalQueue, outWalks);
}

void SGSearchTree::_buildWalksToLeaves(const SGSearchNodePtrDeque& queue, SGWalkVector& outWalks)
{
    for(SGSearchNodePtrDeque::const_iterator iter = queue.begin();
                                             iter != queue.end();
                                             ++iter)
    {
        EdgePtrVec edgeVector;

        // Recursively travel the tree from the leaf to the root collecting the edges in the vector
        addEdgesFromBranch(*iter, edgeVector);

        // Build the walk by adding the edges in reverse order
        SGWalk w(m_pRootNode->getVertex(), false);
        for(EdgePtrVec::reverse_iterator r_iter = edgeVector.rbegin();
                                         r_iter != edgeVector.rend();
                                         ++r_iter)
        {
            w.addEdge(*r_iter);
        }
        outWalks.push_back(w);
    }
}

// Return true if the vertex pX is found somewhere in the branch 
// from pNode to the root
bool SGSearchTree::searchBranchForVertex(SGSearchNode* pNode, Vertex* pX) const
{
    if(pNode == NULL)
        return false;
    if(pNode->getVertex() == pX)
        return true;
    return searchBranchForVertex(pNode->getParent(), pX);
}

//
void SGSearchTree::addEdgesFromBranch(SGSearchNode* pNode, EdgePtrVec& outEdges)
{
    // Terminate the recursion at the root node and dont add an edge
    if(pNode->getParent() != NULL)
    {
        outEdges.push_back(pNode->getEdgeFromParent());
        return addEdgesFromBranch(pNode->getParent(), outEdges);
    }
}

//
void SGSearchTree::_makeFullLeafQueue(SGSearchNodePtrDeque& completeQueue) const
{
    completeQueue.insert(completeQueue.end(), m_expandQueue.begin(), m_expandQueue.end());
    completeQueue.insert(completeQueue.end(), m_goalQueue.begin(), m_goalQueue.end());
    completeQueue.insert(completeQueue.end(), m_doneQueue.begin(), m_doneQueue.end());
}
