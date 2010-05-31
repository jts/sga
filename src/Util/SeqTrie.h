//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SeqTrie.h - Sequence trie data structure
//
#ifndef SEQTRIE_H
#define SEQTRIE_H

#include "Util.h"
#include "Quality.h"
#include <algorithm>
#include <list>

struct PathScore
{
    PathScore() : path_sequence(""), path_corrected(""), 
                  path_score(0.0f), branch_score(0.0f), branch_length(0), branch_cov(0), num_diff(0) {}

    void reverse()
    {
        std::reverse(path_sequence.begin(), path_sequence.end());
        std::reverse(path_corrected.begin(), path_corrected.end());
        std::reverse(probVector.begin(), probVector.end());
    }

    void print()
    {
        printf("CRT: %s PS: %lf BS: %lf AD: %lf\n", path_corrected.c_str(), path_score, branch_score,
                                                    (double)branch_cov / (double)branch_length);
        printf("BSQ: %s\n", path_sequence.c_str());
        printf("QLT: %s\n", Quality::encodeLogProbVector(probVector).c_str());
    }

    // the sequence of nodes that are on this path
    std::string path_sequence;     
    
    // the sequence of the input string, corrected by the path wherever it can be
    std::string path_corrected;     

    double path_score;
    double branch_score;
    int branch_length;
    int branch_cov;
    int num_diff;
    DoubleVector probVector;

};
typedef std::vector<PathScore> PathScoreVector;


class SeqTrie
{
    // Internal datastructures
    class Node;
    struct Link
    {
        // functions
        Link();
        Link(Node* p, char l);
        void increment();
        void decrement();
        void addWeight(double w);

        // data
        Node* pNode;
        char label;
        int count;
        double weight;
    };

    typedef std::list<Link> LinkList;

    class Node
    {
        public:
            // functions
            Node(Node* pParent, char parentLabel);
            ~Node();

            Link* getLink(char label);

            bool insert(const std::string& s, double weight, size_t idx);
            bool remove(const std::string& s, size_t idx);

            void getSequences(std::string curr, StringVector& svOut) const;
            size_t countNodes() const;

            void score(const std::string& s, double lp_correct, 
                                             double lp_error, 
                                             double lp_missing,
                                             size_t idx, const PathScore& curr, PathScoreVector& out);

            void cullChildren(int cutoff);
            void remodel(int cutoff, double weight);

            void writeDot(std::ostream& out) const;

        private:
    
            Link* createChild(char label);
                
            //data
            Link parentLink;
            LinkList childLinks;
    };
    
    // 
    public:

        SeqTrie();
        ~SeqTrie();

        void score(const std::string& s, double p_error, PathScoreVector& out);
        size_t countNodes() const;

        // Creation functions
        void insert(const std::string& s, double weight);
    
        // Remove the string s from the trie
        void remove(const std::string& s);

        // Remove all nodes that have a count less than cutoff from the tree
        void cull(int cutoff);

        // Find links with a count less than cutoff and re-map the branch
        void remodel(int cutoff, double weight);

        // I/O
        void writeDot(std::string filename);

    private:
        
        //
        bool insert(Node* pNode, const std::string& s, size_t idx);
        bool insertAtDepth(Node* pNode, const std::string& s, size_t depth);
        bool remove(Node* pNode, const std::string& s, size_t idx);
        void reap(Node* pNode, int cutoff);

        // Data
        Node* m_pRoot;
};

#endif

