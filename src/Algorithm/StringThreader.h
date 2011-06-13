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
#ifndef STRING_THREADER_H
#define STRING_THREADER_H

#include <list>
#include "BWT.h"
#include "ExtensionDP.h"

// Typedefs
class StringThreaderNode;

typedef std::list<StringThreaderNode*> STNodePtrList;

// A node in the string threading tree
class StringThreaderNode
{

    public:

        // Functions
        StringThreaderNode(const std::string* pQuery, StringThreaderNode* parent);
        ~StringThreaderNode();
        
        // Print the alignment between this branch of the tree and the query
        void printFullAlignment() const;
        
        // 
        void printAllStrings(const std::string& curr) const;

        // Get the alignment error rate over the last context bases
        double getLocalErrorRate(int context) const;
        double getGlobalErrorRate() const;
        int countEdits() const;
        
        // Add a child node to this node with the given label
        // Returns a pointer to the created node
        StringThreaderNode* createChild(const std::string& label);

        // Extend the label of this node by l
        void extend(const std::string& ext);

        // Initialize the alignment columns. This function is for the root node.
        void computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth);
        void computeExtendedAlignment(const std::string& ext, const BandedDPColumn* pPrevColumn);

        //
        bool hasExtensionTerminated() const;

        // Return a suffix of length l of the string represented by this branch
        std::string getSuffix(size_t l) const;

        // Return the complete of this branch including all the parent's labels
        std::string getFullString() const;

    private:
        

        // Data
        // The extension string from the parent
        std::string m_label;

        // The query string being threaded through the graph
        const std::string* m_pQuery;

        // The parent node, can be NULL
        StringThreaderNode* m_pParent;
        STNodePtrList m_children;

        // Alignment information between the label of this node and the query sequence
        // One column per label base
        BandedDPColumnPtrVector m_alignmentColumns;
};

class StringThreader
{
    public:
        StringThreader(const std::string& seed, 
                       const std::string* pQuery,
                       int queryAlignmentEnd,
                       int kmer,
                       const BWT* pBWT, 
                       const BWT* pRevBWT);

        ~StringThreader();

        // Run the threading process
        void run(StringVector& outStrings);
        void printAll();

    private:
        
        // Functions
        void extendLeaves();
        void cullLeavesByLocalError();
        void cullLeavesByEdits();
        
        //
        void checkTerminated(StringVector& outStrings);

        // Perform a 1-base extension of the node
        // Returns true if the node has a branch
        StringVector getDeBruijnExtensions(StringThreaderNode* pNode);

        // Data
        const BWT* m_pBWT; 
        const BWT* m_pRevBWT;
        int m_kmer;
        const std::string* m_pQuery;
        StringThreaderNode* m_pRootNode;
        STNodePtrList m_leaves;
};

#endif
