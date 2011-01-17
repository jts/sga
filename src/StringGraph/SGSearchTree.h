//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SGSearchTree - Data structure used for implementing
// a breadth-first search of a string graph. It is designed
// to return all possible walks between the given start
// and end vertices, up to a given distance.
//
#ifndef SGSEARCHTREE_H
#define SGSEARCHTREE_H

#include "Bigraph.h"
#include "SGWalk.h"
#include <deque>

// foward declare
class SGSearchNode;
typedef std::deque<SGSearchNode*> SGSearchNodePtrDeque;

class SGSearchNode
{
    public:
        SGSearchNode(Vertex* pVertex, EdgeDir expandDir, SGSearchNode* pParent, Edge* pEdgeFromParent);
        ~SGSearchNode();

        // Reduce the child count by 1
        void decrementChildren();

        // Create the children of this node and place pointers to their nodes
        // on the queue. Returns the number of children created;
        int createChildren(SGSearchNodePtrDeque& outQueue);

        SGSearchNode* getParent() const { return m_pParent; }
        Vertex* getVertex() const { return m_pVertex; }
        size_t getDistance() const { return m_extensionDistance; }
        Edge* getEdgeFromParent() const { return m_pEdgeFromParent; }
        int getNumChildren() const { return m_numChildren; }

    private:

        // data
        Vertex* m_pVertex;
        EdgeDir m_expandDir;
        SGSearchNode* m_pParent;
        Edge* m_pEdgeFromParent;

        int m_numChildren;
        size_t m_extensionDistance;
};

class SGSearchTree
{
    public:

        SGSearchTree(Vertex* pStartVertex, 
                     Vertex* pEndVertex,
                     EdgeDir searchDir,
                     size_t distanceLimit,
                     size_t nodeLimit);

        ~SGSearchTree();

        // Returns true if the search has converged on a single vertex. In
        // other words, all walks from the start node share a common vertex,
        // which is represented by one of the nodes waiting expansion.
        // This function is the key to finding collapsed walks that represent
        // complex variation bubbles. It is a heavy operation however, as the
        // entire tree is search for each node in the expand queue. It should
        // only be used on small trees.
        // Returns true if the search converged and the pointer to the vertex
        // is return in pConvergedVertex.
        bool hasSearchConverged(Vertex*& pConvergedVertex);

        // Construct walks representing every path from the start node
        void buildWalksToAllLeaves(SGWalkVector& outWalks);

        // Construct walks representing every path from the start vertex to the goal vertex
        void buildWalksToGoal(SGWalkVector& outWalks);
                
        // Expand all nodes in the queue a single time
        // Returns false if the search has stopped
        bool stepOnce();

        // Check whether the search was aborted or not
        bool wasSearchAborted() const { return m_searchAborted; }

    private:

        // Search the branch from pNode to the root for pX.  
        bool searchBranchForVertex(SGSearchNode* pNode, Vertex* pX) const;

        // Build the walks from the root to the leaves in the queue
        void _buildWalksToLeaves(const SGSearchNodePtrDeque& queue, SGWalkVector& outWalks);
        void addEdgesFromBranch(SGSearchNode* pNode, EdgePtrVec& outEdges);

        // Build a queue with all the leaves in it
        void _makeFullLeafQueue(SGSearchNodePtrDeque& completeQueue) const;

        // We keep the pointers to the search nodes
        // in one of three queues. 
        // The goal queue contains the nodes representing the vertex we are searching for.
        // The expand queue contains nodes that have not yet been explored.
        // The done queue contains non-goal nodes that will not be expanded further.
        // Together, they represent all leaves of the tree
        SGSearchNodePtrDeque m_goalQueue;
        SGSearchNodePtrDeque m_expandQueue;
        SGSearchNodePtrDeque m_doneQueue;
        size_t m_totalNodes; // The total number of nodes in the search tree
    
        SGSearchNode* m_pRootNode;
        Vertex* m_pGoalVertex;
        EdgeDir m_initialDir;

        size_t m_distanceLimit;
        size_t m_nodeLimit;

        // Flag indicating the search was aborted
        bool m_searchAborted;

};

#endif
