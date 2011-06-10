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

// Return a suffix of length l of the string represented by this branch
std::string StringThreaderNode::getFullString() const
{
    if(pParent == NULL)
        return label;
    else
        return pParent->getFullString() + label;
}

// Create a child node
StringThreaderNode* StringThreaderNode::createChild(const std::string& label)
{
    StringThreaderNode* pAdded = new StringThreaderNode(label, this);
    m_children.push_back(pAdded);
    return pAdded;
}

// Extend the label of this node
void StringThreaderNode::extend(const std::string& ext)
{
    label.append(ext);
}

//
void StringThreaderNode::printFullAlignment(const std::string* pQuery) const
{
    std::string fullString = getFullString();
    StdAlnTools::printGlobalAlignment(fullString, *pQuery);
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

    std::string test1 = "CGTGCGATGCACTGCTACGGCTCGCCTAGATCA";
    std::string test2 = "CGTGCGATGCACTGCATACGGCTCGCCTAGATCA";
    GlobalAlnParams params;
//    ExtensionDP::initialize(seed, pQuery->substr(0, seed.size()), params);
    ExtensionDP::initialize(test1, test2, params);

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
    int count = 100;
    while(count-- > 0)
        extendLeaves();
}

// Extend each leaf node
void StringThreader::extendLeaves()
{
    std::cout << "****Extension***\n\n";
    STNodePtrList newLeaves;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        std::cout << "LEAF BEFORE EXTENSION\n";
        (*iter)->printFullAlignment(m_pQuery);

        StringVector extensions = getDeBruijnExtensions(*iter);

        // Either extend the current node or branch it
        // If no extension, do nothing so this node
        // is no longer marked a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front());
            newLeaves.push_back(*iter);
        }
        else if(extensions.size() > 1)
        {
            // Branch
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                StringThreaderNode* pAdded = (*iter)->createChild(extensions[i]);
                newLeaves.push_back(pAdded);
            }
        }
    }

    m_leaves = newLeaves;
}

// Extend a leaf node, possibly creating a branch
StringVector StringThreader::getDeBruijnExtensions(StringThreaderNode* pNode)
{
    // Get the last k-1 bases of the node
    std::string pmer = pNode->getSuffix(m_kmer - 1);
    assert((int)pmer.size() == m_kmer - 1);
    std::string rc_pmer = reverseComplement(pmer);

    // Get an interval for the p-mer and its reverse complement
    BWTIntervalPair ip = BWTAlgorithms::findIntervalPair(m_pBWT, m_pRevBWT, pmer);
    BWTIntervalPair rc_ip = BWTAlgorithms::findIntervalPair(m_pBWT, m_pRevBWT, rc_pmer);

    // Get the extension bases
    AlphaCount64 extensions;
    AlphaCount64 rc_extensions;
    if(ip.interval[1].isValid())
        extensions += BWTAlgorithms::getExtCount(ip.interval[1], m_pRevBWT);
    if(rc_ip.interval[1].isValid())
        rc_extensions = BWTAlgorithms::getExtCount(rc_ip.interval[0], m_pBWT);
    rc_extensions.complement();
    extensions += rc_extensions;

    // Loop over the DNA symbols, if there is are more than two characters create a branch
    // otherwise just perform an extension.
    bool hasExtension = extensions.hasDNAChar();

    StringVector out;
    if(hasExtension)
    {
        for(int i = 0; i < DNA_ALPHABET::size; ++i)
        {
            char b = DNA_ALPHABET::getBase(i);
            if(extensions.get(b) > 0)
            {
                // extend to b
                std::string tmp;
                tmp.append(1,b);
                out.push_back(tmp);
            }
        }
    }

    // If the node branched, return true so the outer function can remove it from the leaf list
    return out;
}
