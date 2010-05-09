#ifndef SUFFIXTREE_H
#define SUFFIXTREE_H

#include <assert.h>
#include <fstream>
#include "STCommon.h"


//
// Typedefs 
//
typedef std::string STLabel;
typedef uint8_t STIdx;

const uint8_t MAX_CHILDREN = 5;

struct STNode
{
    STNode(STLabel l, STNode* p) : edgeLabel(l)
    {
        parent = p;
        for(STIdx i = 0; i < MAX_CHILDREN; ++i)
        {
            children[i] = NULL;
        }
    }

    bool isLeaf() const
    {
        if(edgeLabel.size() > 0)
            return edgeLabel[edgeLabel.length() - 1] == '$';
        else
            return false;
    }

    // Returns the (approximate?) size of the node
    size_t getByteSize() const
    {
        size_t nodeSize = sizeof(STNode);
        size_t strSize = edgeLabel.size();
        //std::cout << "NS\t" << nodeSize << "\tSS\t" << strSize << "\n";
        return nodeSize + strSize;
    }

    STLabel edgeLabel; // The label from the parent to this node
    STNode* parent;
    STNode* children[MAX_CHILDREN];
};


struct PathEnd
{
    PathEnd(STNode* p, size_t em, size_t po) : pNode(p), edgeMatch(em), prefixOffset(po) {}
    
    STNode* pNode; // The node a path ended on
    size_t edgeMatch; // The number of characters along its edge that were matched
    size_t prefixOffset; // The position of the prefix that the match along the edge starts
};

struct STLeaf
{
    int index;
};

class SuffixTree
{
    public:
        SuffixTree(std::string str);
        ~SuffixTree();
        // Insert a string into the suffix tree
        void insert(std::string s);

        // Write the suffix tree as a dot file
        void printInfo() const;
        void writeDot(std::string filename) const;
        size_t getByteSize() const;
    
    private:
        
        // Insert a single suffix into the tree
        void insertSuffix(std::string s);

        // Print suffixes
        void printSuffixes(STNode* pNode) const;

        // Print a single suffix, beginning for the given node
        void printSuffix(STNode* pNode) const;
        
        // Count the number of nodes in the subtree
        size_t count(STNode* pNode) const;

        // Construct the tree from string s
        void construct(std::string s);

        // Find the path from the root that maximally matches a prefix of s
        PathEnd findPath(std::string s) const;

        // Return the longest matching prefix of s,t
        size_t prefixMatch(std::string s, std::string t) const;

        // Validate the pointers in the tree
        void validate(STNode* pNode) const;

        // Write the node in dot format
        void writeNode(std::ostream& out, STNode* pNode) const;

        // Count the size of the node and its children
        size_t getByteSize(STNode* pNode) const;

        // Data
        STNode* m_pRoot;
        size_t m_suffixesInserted;
};

#endif

