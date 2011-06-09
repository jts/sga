///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringThreader - Iteratively construct a
// set of strings by threading a query sequence
// through a graph.
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "StringThreader.h"
#include "BWTAlgorithms.h"
#include "LRAlignment.h"
#include "StdAlnTools.h"

//
// StringThreaderNode
//
StringThreaderNode::StringThreaderNode(const std::string& l, StringThreaderNode* parent) : label(l),
                                                                                           pParent(parent)

{


}

// Delete the children of the node
StringThreaderNode::~StringThreaderNode()
{
    for(STNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;
}

// Return a suffix of length l of the string represented by this branch
std::string StringThreaderNode::getSuffix(size_t l) const
{
    size_t n = label.size();
    if(l <= n)
    {
        return label.substr(n - l, l);
    }
    else 
    {
        assert(pParent != NULL);
        return pParent->getSuffix(l - n) + label;
    }
}

// Create a child node
StringThreaderNode* StringThreaderNode::createChild(const std::string& label)
{
    StringThreaderNode* pAdded = new StringThreaderNode(label, this);
    m_children.push_back(pAdded);
    return pAdded;
}

//
void StringThreaderNode::printFullAlignment(const std::string* pQuery) const
{
    StdAlnTools::printGlobalAlignment(label, *pQuery);
}

//
// StringTheader
//
StringThreader::StringThreader(const std::string& seed, 
                               const std::string* pQuery, 
                               int kmer, 
                               const BWT* pBWT, 
                               const BWT* pRevBWT) : m_pBWT(pBWT), m_pRevBWT(pRevBWT), 
                                                     m_kmer(kmer), m_pQuery(pQuery)
{
    // Create the root node containing the seed string
    m_pRootNode = new StringThreaderNode(seed, NULL);
    m_leaves.push_back(m_pRootNode);
}

//
StringThreader::~StringThreader()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

// Run the threading algorithm
void StringThreader::run()
{
    // Extend the leaf nodes
    extendLeaves();
}

// Extend each leaf node
void StringThreader::extendLeaves()
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        performExtension(*iter);
    }
}

// Extend a leaf node, possibly creating a branch
void StringThreader::performExtension(StringThreaderNode* pNode)
{
    // Get the last k-1 bases of the node
    std::string pmer = pNode->getSuffix(m_kmer - 1);
    assert(!pmer.empty());
    std::cout << "PMER: " << pmer << "\n";
}
