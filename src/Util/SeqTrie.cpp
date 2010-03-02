//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SeqTrie.h - Sequence trie data structure
//
#include "SeqTrie.h"
#include <iostream>
#include <fstream>

// Link
void SeqTrie::Link::increment()
{
	++count;
}

void SeqTrie::Link::decrement()
{
	--count;
}

// Node
SeqTrie::Node::Node(Node* pParent, char parentLabel)
{
	parentLink.pNode = pParent;
	parentLink.label = parentLabel;
}

// Destructor, destroy the children of this node
SeqTrie::Node::~Node()
{
	for(LinkList::iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
		delete iter->pNode;
}

// Return a pointer to the link with label otherwise NULL
SeqTrie::Link* SeqTrie::Node::getLink(char label)
{
	for(size_t i = 0; i < pChildLinks.size(); ++i)
	{
		for(LinkList::iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
			if(iter->label == label)
				return &(*iter);
	}
	return NULL;
}

// Create a new child node
SeqTrie::Link* SeqTrie::Node::createChild(char label)
{
	Node* pChild = new Node(this, label);
	Link l = {pChild, label, 0};
	pChildLinks.push_back(l);
	return &pChildLinks.back();
}

// Remove children with count below cutoff
void SeqTrie::Node::removeChildren(int cutoff)
{
	LinkList::iterator iter = pChildLinks.begin(); 
	while(iter != pChildLinks.end())
	{
		if(iter->count < cutoff)
		{
			delete iter->pNode; // recursive
			iter = pChildLinks.erase(iter);
		}
		else
		{
			iter->pNode->removeChildren(cutoff);
			++iter;
		}
	}
}

// Recursive dot writer function
void SeqTrie::Node::writeDot(std::ostream& out) const
{
	out << "\"" << this << "\" [label=\"\"];\n";
	for(LinkList::const_iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
	{
		out << "\"" << this << "\" -> \"" << iter->pNode << "\" [label=\""
		    << iter->label << "," << iter->count << "\"];\n";
		iter->pNode->writeDot(out);
	}
}

//
SeqTrie::SeqTrie()
{
	m_pRoot = new Node(NULL, '^');
}

// insert the string s into the trie
void SeqTrie::insert(const std::string& s)
{
	if(s.empty())
		return;
	insert(m_pRoot, s, 0);
}

// remove s from the trie
void SeqTrie::remove(const std::string& s)
{
	WARN_ONCE("SeqTrie::remove only decrementing but not reaping");
	remove(m_pRoot, s, 0);
}

// remove all sub-tries that have a link less than cutoff
void SeqTrie::reap(int cutoff)
{
	m_pRoot->removeChildren(cutoff);
}

// insert the string s into the trie so that it is a child 
// of the node(s) at DEPTH. Children of the root (the first
// nodes) are depth 0. If depth is higher than the deepest
// node in the trie, this will do nothing.
void SeqTrie::insertAtDepth(const std::string& s, size_t depth)
{
	bool inserted = insertAtDepth(m_pRoot, s, depth);
	if(!inserted)
	{
		std::cerr << "Error: Could not insert string at depth " << depth << "\n";
		assert(false);
	}
}

// Recursive insertDepth call, descend into child branches into the depth
// insertion point has been found then call insert
bool SeqTrie::insertAtDepth(Node* pNode, const std::string& s, size_t depth)
{
	if(depth == 0)
	{
		return insert(pNode, s, 0);
	}
	else
	{
		bool rc = false;
		for(LinkList::iterator iter = pNode->pChildLinks.begin(); iter != pNode->pChildLinks.end(); ++iter)
			rc = insertAtDepth(iter->pNode, s, depth - 1) || rc;
		return rc;
	}
}

// insert the string s into the trie starting at pNode and character idx
// returns true on successfull insert
bool SeqTrie::insert(Node* pNode, const std::string& s, size_t idx)
{
	assert(pNode != NULL);
	char b = s[idx];
	Link* pLink = pNode->getLink(b);
	if(pLink == NULL)
	{
		// The node does not have a link with this label, create a new node
		pLink = pNode->createChild(b);
	}
	pLink->increment();

	// Recurse
	if(++idx != s.size())
		return insert(pLink->pNode, s, idx);
	else
		return true;
}

// remove the string s from the trie starting at pNode and character idx
bool SeqTrie::remove(Node* pNode, const std::string& s, size_t idx)
{
	assert(pNode != NULL);
	char b = s[idx];
	Link* pLink = pNode->getLink(b);
	if(pLink != NULL)
	{
		pLink->decrement();
		if(++idx != s.size())
			return remove(pLink->pNode, s, idx);
	}
	return false;
}

// Write the trie to a dot file
void SeqTrie::writeDot(std::string filename)
{
	std::ofstream writer(filename.c_str());
	writer << "digraph G\n{\n";
	m_pRoot->writeDot(writer);
	writer << "}\n";
	writer.close();
}

//
SeqTrie::~SeqTrie()
{
	delete m_pRoot; // recursively destroys children
	m_pRoot = NULL;
}
