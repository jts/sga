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

STIdx base2Idx(char b);

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
		void writeDot(std::string filename);
	
	private:
		
		// Count the number of nodes in the subtree
		size_t count(STNode* pNode);

		// Construct the tree from string s
		void construct(std::string s);

		// Find the path from the root that maximally matches a prefix of s
		PathEnd findPath(std::string s);

		// Return the longest matching prefix of s,t
		size_t prefixMatch(std::string s, std::string t);

		// Validate the pointers in the tree
		void validate(STNode* pNode);

		// Write the node in dot format
		void writeNode(std::ostream& out, STNode* pNode);

		STNode* m_pRoot;
};

#endif

