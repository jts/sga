#include "STCommon.h"
#include "SuffixTree.h"

//
//
//
SuffixTree::SuffixTree(std::string str)
{
	construct(str);
	validate(m_pRoot);
}

//
// Naive O(n^2) suffix tree construction algorithm from Gusfield
//
void SuffixTree::construct(std::string s) 
{
	// Create the new root and insert the first suffix, consisting of the entire string
	m_pRoot = new STNode("", NULL);
	m_pRoot->children[getBaseRank(s[0])] = new STNode(s, m_pRoot);
	m_suffixesInserted = 1;

	// Insert the rest of the suffixes
	insert(s.substr(1));
	printInfo();
	writeDot("st.dot");
	printSuffixes(m_pRoot);
}

//
// Insert a string into the tree
//
void SuffixTree::insert(std::string s)
{
	for(size_t i = 0; i < s.length(); ++i)
		insertSuffix(s.substr(i));
	printSuffixes(m_pRoot);
}

//
// Print information about the tree
//
void SuffixTree::printInfo() const
{
	size_t l = m_suffixesInserted;
	size_t numNodes = count(m_pRoot);
	size_t b = getByteSize();

	std::cout << "Suffix tree -- Total string length: " << l
			<< " nodes: " << numNodes
			<< " bytes: " << b 
			<< " bpc: " << (double)b / (double)l << "\n";
}

//
// insert a suffix into the suffix tree
//
void SuffixTree::insertSuffix(std::string s)
{
	// Match characters from S[i..len] along the tree
	PathEnd pe = findPath(s);
	std::cout << "Substring " << s << " ends in node " << pe.pNode << " em: " << pe.edgeMatch << " po: " << pe.prefixOffset << "\n";

	std::string matchedStr = s.substr(pe.prefixOffset, pe.edgeMatch);
	std::string unmatchedStr = s.substr(pe.prefixOffset + pe.edgeMatch);

	// If the unmatched length of the string is 0 it means the suffix is already in the tree
	// Nothing needs to be done
	if(unmatchedStr.length() == 0)
	{
		return;
	}

	// Case 1: The match does not break an edge
	if(pe.edgeMatch == 0)
	{
		// Create a new node, from the current node with the unmatched part of curr as the label
		pe.pNode->children[getBaseRank(unmatchedStr[0])] = new STNode(unmatchedStr, pe.pNode);
	}
	else // Case 2: The match breaks an edge
	{
		// Split the edge into the matched and unmatched parts
		std::string matchedEdge = pe.pNode->edgeLabel.substr(0, pe.edgeMatch);
		std::string unmatchedEdge = pe.pNode->edgeLabel.substr(pe.edgeMatch);

		// Create the new node and add it to the current nodes parent
		std::cout << "MatchedCurr:\t" << matchedStr << "\n";
		std::cout << "UnmatchedStr:\t" << unmatchedStr << "\n";
		std::cout << "MatchedEdge:\t" << matchedEdge << "\n";
		std::cout << "UnmatchedEdge:\t" << unmatchedEdge << "\n";
		
		// Create a new node, consisting of the matched portion of curr and the edge
		// Insert it between pe.pNode and the parent
		STNode* pBridge = new STNode(matchedStr, pe.pNode->parent);
		assert(matchedStr == matchedEdge);
		pe.pNode->parent->children[ getBaseRank(matchedStr[0]) ] = pBridge;

		// Add the unmatched portion the edge to the newly created node
		assert(unmatchedEdge[0] != unmatchedStr[0]); // if this was true we could extend the match further
		pe.pNode->edgeLabel = unmatchedEdge; // update the label of the split node
		pe.pNode->parent = pBridge;
		pBridge->children[ getBaseRank(unmatchedEdge[0]) ] = pe.pNode; // add it to the bridge node

		// Create a node consisting of the unmatched portion of this suffix and add it
		pBridge->children[ getBaseRank(unmatchedStr[0]) ] = new STNode(unmatchedStr, pBridge);
	}
	m_suffixesInserted++;
}

//
// Return the size of the suffix tree in bytes
//
size_t SuffixTree::getByteSize() const
{
	return getByteSize(m_pRoot);
}

//
// Find the path from the root of the tree that matches a prefix in s
// Return the node
//
PathEnd SuffixTree::findPath(std::string s) const
{
	size_t len = s.size();
	STNode* pCurr = m_pRoot;
	bool done = false;
	size_t i = 0;
	while(!done && i < len)
	{
		STIdx childIdx = getBaseRank(s[i]);
		
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
size_t SuffixTree::prefixMatch(std::string s, std::string t) const
{
	size_t stop = (s.length() < t.length()) ? s.length() : t.length();
	size_t i = 0;
	while(s[i] == t[i] && i < stop)
		++i;
	return i;
}

//
// Print all the suffixes held in the tree in nodes below pNode
//
void SuffixTree::printSuffixes(STNode* pNode) const
{
	if(pNode->isLeaf())
	{
		printSuffix(pNode);
	}
	else
	{
		for(STIdx i = 0; i < MAX_CHILDREN; ++i)
		{
			if(pNode->children[i] != NULL)
				printSuffixes(pNode->children[i]);
		}
	}
}

//
// Print a suffix, starting from the given node
//
void SuffixTree::printSuffix(STNode* pNode) const
{
	STLabel str;
	while(pNode != m_pRoot)
	{
		str = pNode->edgeLabel + str;
		pNode = pNode->parent;
	}
	std::cout << str << "\n";
}

//
// Count the number of nodes in the subtree
//
size_t SuffixTree::count(STNode* pNode) const
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
// Get the size of the subtree (in bytes) 
// 
size_t SuffixTree::getByteSize(STNode* pNode) const
{
	size_t size = 0;
	for(STIdx i = 0; i < MAX_CHILDREN; ++i)
	{
		if(pNode->children[i] != NULL)
			size += getByteSize(pNode->children[i]);
	}
	return size + pNode->getByteSize();
}


//
// Validate the links in the tree
//
void SuffixTree::validate(STNode* pNode) const
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
void SuffixTree::writeDot(std::string filename) const
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
void SuffixTree::writeNode(std::ostream& out, STNode* pNode) const
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
