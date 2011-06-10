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
StringThreaderNode::StringThreaderNode(const std::string* pQuery,
                                       StringThreaderNode* parent) : m_pQuery(pQuery),
                                                                     m_pParent(parent)

{

}

// Delete the children of the node
StringThreaderNode::~StringThreaderNode()
{
    // Delete children
    for(STNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;

    // Delete alignment columns
    for(size_t i = 0; i < m_alignmentColumns.size(); ++i)
        delete m_alignmentColumns[i];
}

// Return a suffix of length l of the string represented by this branch
std::string StringThreaderNode::getSuffix(size_t l) const
{
    size_t n = m_label.size();
    if(l <= n)
    {
        return m_label.substr(n - l, l);
    }
    else 
    {
        assert(m_pParent != NULL);
        return m_pParent->getSuffix(l - n) + m_label;
    }
}

// Return a suffix of length l of the string represented by this branch
std::string StringThreaderNode::getFullString() const
{
    if(m_pParent == NULL)
        return m_label;
    else
        return m_pParent->getFullString() + m_label;
}

// Create a child node
StringThreaderNode* StringThreaderNode::createChild(const std::string& label)
{
    StringThreaderNode* pAdded = new StringThreaderNode(m_pQuery, this);
    assert(!m_alignmentColumns.empty());
    pAdded->computeExtendedAlignment(label, m_alignmentColumns.back());
    m_children.push_back(pAdded);
    return pAdded;
}

// Extend the label of this node
void StringThreaderNode::extend(const std::string& ext)
{
    assert(!ext.empty());
    assert(!m_alignmentColumns.empty());
    computeExtendedAlignment(ext, m_alignmentColumns.back());
}
   
void StringThreaderNode::computeExtendedAlignment(const std::string& ext, const BandedDPColumn* pPrevColumn)
{
    // Calculate a new alignment column for each base of the label
    for(size_t i = 0; i < ext.size(); ++i)
    {
        BandedDPColumn* pNewColumn = ExtensionDP::createNewColumn(ext[i], *m_pQuery, pPrevColumn);
        m_alignmentColumns.push_back(pNewColumn);
        pPrevColumn = pNewColumn;
    }
    m_label.append(ext);
}

// Initialize the alignment columns
void StringThreaderNode::computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth)
{
    // Create the initial alignment columnds between label and query
    m_label = initialLabel;
    assert(!m_label.empty());
    ExtensionDP::createInitialAlignment(m_label, m_pQuery->substr(0, queryAlignmentEnd), bandwidth, m_alignmentColumns);
}

//
double StringThreaderNode::getLocalErrorRate(int context) const
{
    return ExtensionDP::calculateLocalEditPercentage(m_alignmentColumns.back(), context);
}

double StringThreaderNode::getGlobalErrorRate() const
{
    return ExtensionDP::calculateGlobalEditPercentage(m_alignmentColumns.back());
}

//
void StringThreaderNode::printFullAlignment() const
{
    std::string fullString = getFullString();

//    std::cout << "STDALN ALIGNMENT:\n";
//    StdAlnTools::printGlobalAlignment(*m_pQuery, fullString);

//    std::cout << "EXTENSIONDP ALIGNMENT:\n";
    ExtensionDP::printAlignment(fullString, *m_pQuery, m_alignmentColumns.back());

    double localER = ExtensionDP::calculateLocalEditPercentage(m_alignmentColumns.back(), 20);
    printf("LocalER: %lf\n", localER);
}

//
void StringThreaderNode::printAllStrings(const std::string& curr) const
{
    if(m_children.empty())
    {
        std::cout << "S: " << curr + m_label << "\n";
    }
    else
    {
        for(STNodePtrList::const_iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
            (*iter)->printAllStrings(curr + m_label);
    }
}

//
// StringTheader
//
StringThreader::StringThreader(const std::string& seed, 
                               const std::string* pQuery,
                               int queryAlignmentEnd,
                               int kmer, 
                               const BWT* pBWT, 
                               const BWT* pRevBWT) : m_pBWT(pBWT), m_pRevBWT(pRevBWT), 
                                                     m_kmer(kmer), m_pQuery(pQuery)
{
    // Create the root node containing the seed string
    m_pRootNode = new StringThreaderNode(pQuery, NULL);
    m_pRootNode->computeInitialAlignment(seed, queryAlignmentEnd, 50);
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
    while(!m_leaves.empty())
        extendLeaves();
    printAll();
}

// Print the string represented by every node
void StringThreader::printAll()
{
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node
void StringThreader::extendLeaves()
{
    STNodePtrList newLeaves;
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
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

    m_leaves.clear();

    int context = 20;
    double threshold = 0.3f;
    // Calculate the local error rate of the alignments to each new leaf
    // If it is less than threshold, add the leaf to the node
    std::cout << "***EXTENSION:\n";
    for(STNodePtrList::iterator iter = newLeaves.begin(); iter != newLeaves.end(); ++iter)
    {
        double ler = (*iter)->getLocalErrorRate(context);
        double ger = (*iter)->getGlobalErrorRate();
        printf("ger: %.2lf\n", ger);

        if(ler < threshold)
        {
            printf("Keeping leaf with ler: %.2lf ger: %.2lf\n", ler, ger);
            m_leaves.push_back(*iter);
        }
        else
        {
            printf("Culling leaf with ler: %.2lf ger: %.2lf\n", ler, ger);
            
            /*
            static int boom = 20;
            if(boom-- == 0)
                exit(1);
            */
        }
        (*iter)->printFullAlignment();
    }

    std::cout << "Leaves now: " << m_leaves.size() << "\n";
    if(m_leaves.size() > 100)
        exit(1);
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
