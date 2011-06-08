///-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// StringBuilder - Iteratively construct strings
// that represent sequences in an assembly graph.
// The assembly graph is abstractly represented as
// an FM-index.
//
#ifndef STRING_BUILDER_H
#define STRING_BUILDER_H

#include <list>
#include "BWT.h"

// A node in the string threading tree
class StringThreaderNode
{
    public:

        // Functions
        StringThreaderNode(const std::string& l, int tae, int qae, int as);
        ~StringThreaderNode();

        void printFullAlignment(const std::string* pQuery) const;

        // Data
        // The extension string from the parent
        std::string label;

        // The parent node, can be NULL
        StringThreaderNode* pParent;
        typedef std::list<StringThreaderNode*> STNodePtrList;
        STNodePtrList m_children;

        // The coordinates of the end of the alignment
        // on the query string and the target thread up to this point
        int target_alignment_end;
        int query_alignment_end;

        // The score of the alignment string
        int alignment_score;
};

class StringThreader
{
    public:
        StringThreader(const std::string& seed, 
                       const std::string* pQuery, 
                       int kmer,
                       const BWT* pBWT, 
                       const BWT* pRevBWT);
        ~StringThreader();

        void run();

    private:
        
        const BWT* m_pBWT; 
        const BWT* m_pRevBWT;
        int m_kmer;
        const std::string* m_pQuery;
        StringThreaderNode* m_pRootNode;
};

class StringBuilder
{
    public:
        StringBuilder(const std::string& seed, int kmer, const BWT* pBWT, const BWT* pRevBWT);

        // Perform one round of extension for all the strings
        void extendOnce();

        // Remove all strings but the top n scoring
        void cull(const std::string& query, size_t n);

        // Print the set of strings created
        void print() const;
        void print(const std::string& query) const;

    private:

        //
        const BWT* m_pBWT; 
        const BWT* m_pRevBWT;

        StringVector m_strings;
        int m_kmer;
};

#endif
