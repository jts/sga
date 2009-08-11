#include "SuffixTree.h"

STIdx base2Idx(char b)
{
	switch(b)
	{
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		case '$':
			return 4;
		default:
			assert(false);
	}
}

//
//
//
SuffixTree::SuffixTree(std::string str)
{
	construct(str);
	validate(m_pRoot);
	size_t numNodes = count(m_pRoot);
	std::cout << "Constructed a tree from a string of length " << str.length() << " with " << numNodes << " nodes\n";
}

//
// Naive O(n^2) suffix tree construction algorithm from Gusfield
//
void SuffixTree::construct(std::string s)
{
	// Create the new root and insert the first suffix, consisting of the entire string
	m_pRoot = new STNode("", NULL);
	m_pRoot->children[base2Idx(s[0])] = new STNode(s, m_pRoot);

	size_t len = s.length();
	std::cout << "ROOT: " << m_pRoot << "\n";
	for(size_t i = 1; i < len; ++i)
	{
		std::string curr = s.substr(i);
		// Match characters from S[i..len] along the tree
		PathEnd pe = findPath(curr);
		std::cout << "Substring " << curr << " ends in node " << pe.pNode << " em: " << pe.edgeMatch << " po: " << pe.prefixOffset << "\n";

		std::string matchedCurr = curr.substr(pe.prefixOffset, pe.edgeMatch);
		std::string unmatchedCurr = curr.substr(pe.prefixOffset + pe.edgeMatch);
		assert(unmatchedCurr.length() > 0); // At the very least the $ must be unmatched

		// Case 1: The match does not break an edge
		if(pe.edgeMatch == 0)
		{
			// Create a new node, from the current node with the unmatched part of curr as the label
			pe.pNode->children[base2Idx(unmatchedCurr[0])] = new STNode(unmatchedCurr, pe.pNode);
		}
		else // Case 2: The match breaks an edge
		{
			// Split the edge into the matched and unmatched parts
			std::string matchedEdge = pe.pNode->edgeLabel.substr(0, pe.edgeMatch);
			std::string unmatchedEdge = pe.pNode->edgeLabel.substr(pe.edgeMatch);

			// Create the new node and add it to the current nodes parent
			std::cout << "MatchedCurr:\t" << matchedCurr << "\n";
			std::cout << "UnmatchedCurr:\t" << unmatchedCurr << "\n";
			std::cout << "MatchedEdge:\t" << matchedEdge << "\n";
			std::cout << "UnmatchedEdge:\t" << unmatchedEdge << "\n";
			
			// Create a new node, consisting of the matched portion of curr and the edge
			// Insert it between pe.pNode and the parent
			STNode* pBridge = new STNode(matchedCurr, pe.pNode->parent);
			assert(matchedCurr == matchedEdge);
			pe.pNode->parent->children[ base2Idx(matchedCurr[0]) ] = pBridge;

			// Add the unmatched portion the edge to the newly created node
			assert(unmatchedEdge[0] != unmatchedCurr[0]); // if this was true we could extend the match further
			pe.pNode->edgeLabel = unmatchedEdge; // update the label of the split node
			pe.pNode->parent = pBridge;
			pBridge->children[ base2Idx(unmatchedEdge[0]) ] = pe.pNode; // add it to the bridge node

			// Create a node consisting of the unmatched portion of this suffix and add it
			pBridge->children[ base2Idx(unmatchedCurr[0]) ] = new STNode(unmatchedCurr, pBridge);
		}

		writeDot("st.dot");
	}
}

//
// Find the path from the root of the tree that matches a prefix in s
// Return the node
//
PathEnd SuffixTree::findPath(std::string s)
{
	size_t len = s.size();
	STNode* pCurr = m_pRoot;
	bool done = false;
	size_t i = 0;
	while(!done && i < len)
	{
		STIdx childIdx = base2Idx(s[i]);
		
		// If the current node does not have a child starting with this letter
		if(pCurr->children[childIdx] != NULL)
		{
			// Match from [i..l] along the edge of the child
			const std::string& elabel = pCurr->children[childIdx]->edgeLabel;
			size_t matchSize = prefixMatch(s.substr(i), elabel);

			std::cout << "Matchsize: " << s.substr(i) << " " << elabel << " = " << matchSize << "\n";

			// Check if the entire prefix has now been matched or if the edge is split
			if(i + matchSize == len || matchSize != elabel.size())
			{
				return PathEnd(pCurr->children[childIdx], matchSize, i);
			}
			else
			{
				pCurr = pCurr->children[childIdx];
				i += matchSize;
			}
		}
		else
		{
			std::cout << "Node does not have a path starting with " << s[i] << "\n";
			return PathEnd(pCurr, 0, i);
		}
	}
	assert(false);
	return PathEnd(NULL, 0, i);
}

//
// Match the prefixes of two strings returning the number of characters matched
//
size_t SuffixTree::prefixMatch(std::string s, std::string t)
{
	size_t s_len = s.length();
	size_t t_len = t.length();
	size_t stop = (s_len < t_len) ? s_len : t_len;
	size_t i = 0;
	while(s[i] == t[i] && i < stop)
		++i;
	return i;
}

//
// Count the number of nodes in the subtree
//
size_t SuffixTree::count(STNode* pNode)
{
	size_t childCount = 0;
	for(STIdx i = 0; i < MAX_CHILDREN; ++i)
	{
		if(pNode->children[i] != NULL)
			childCount += count(pNode->children[i]);
	}
	return childCount + 1;
}

//
// Validate the links in the tree
//
void SuffixTree::validate(STNode* pNode)
{
	for(STIdx i = 0; i < MAX_CHILDREN; ++i)
	{
		if(pNode->children[i] != NULL)
		{
			assert(pNode->children[i]->parent == pNode);
			validate(pNode->children[i]);
		}
	}
}

//
//
//
void SuffixTree::writeDot(std::string filename)
{
	std::ofstream out(filename.c_str());
	
	out << "digraph G\n{\n";

	writeNode(out, m_pRoot);

	out << "}\n";
	out.close();
}

//
//
//
void SuffixTree::writeNode(std::ostream& out, STNode* pNode)
{
	out << "\"" << pNode << "\" [ label =\"" << pNode << "\"]\n";
	for(STIdx i = 0; i < MAX_CHILDREN; ++i)
	{
		if(pNode->children[i] != NULL)
		{
			out << "\"" << pNode << "\" -> \"" << pNode->children[i] << "\" [ label=\"" << pNode->children[i]->edgeLabel << "\"]\n";
			writeNode(out, pNode->children[i]);
		}
	}
}
