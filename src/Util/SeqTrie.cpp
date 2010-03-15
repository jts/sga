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
#include <math.h>

//
// Link
//
SeqTrie::Link::Link() : pNode(NULL), label('\0'), count(0), weight(0.0f)
{

}

SeqTrie::Link::Link(Node* p, char l) : pNode(p), label(l), count(0), weight(0.0f)
{

}

void SeqTrie::Link::increment()
{
	++count;
}

void SeqTrie::Link::decrement()
{
	--count;
}

void SeqTrie::Link::addWeight(double w)
{
	weight += w;
}

//
// Node
//
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
	for(LinkList::iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
		if(iter->label == label)
			return &(*iter);
	return NULL;
}

// Create a new child node
SeqTrie::Link* SeqTrie::Node::createChild(char label)
{
	Node* pChild = new Node(this, label);
	Link l(pChild, label);
	pChildLinks.push_back(l);
	return &pChildLinks.back();
}

// Score the string s against the trie, descending from this node
void SeqTrie::Node::score(const std::string& s, double p_error, size_t idx, const PathScore& curr, PathScoreVector& out)
{
	if(s.size() == idx || pChildLinks.empty())
	{
		out.push_back(curr);
	}

	// Descend into each subtree
	for(LinkList::iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
	{
		PathScore branch = curr;
		if(iter->label == s[idx])
			branch.score += log(1 - p_error);
		else
			branch.score += log(p_error);
		branch.path.append(1, iter->label);
		iter->pNode->score(s, p_error, idx + 1, branch, out);
	}
}


// insert the string s into this node starting with the symbol at idx
bool SeqTrie::Node::insert(const std::string& s, double weight, size_t idx)
{
	char b = s[idx];
	Link* pLink = getLink(b);
	if(pLink == NULL)
	{
		pLink = createChild(b);
	}
	
	pLink->increment();
	pLink->addWeight(weight);

	// Recurse
	if(++idx != s.size())
		return pLink->pNode->insert(s, weight, idx);
	else
		return true;
}

// remove the string s from the trie starting at pNode and character idx
bool SeqTrie::Node::remove(const std::string& s, size_t idx)
{
	char b = s[idx];
	Link* pLink = getLink(b);
	if(pLink != NULL)
	{
		pLink->decrement();
		if(++idx != s.size())
			return pLink->pNode->remove(s, idx);
	}
	return false;
}


// Remove children with count below cutoff
void SeqTrie::Node::cullChildren(int cutoff)
{
	LinkList::iterator iter = pChildLinks.begin(); 
	while(iter != pChildLinks.end())
	{
		if(pChildLinks.size() > 1 && iter->count < cutoff)
		{
			delete iter->pNode; // recursive
			iter = pChildLinks.erase(iter);
		}
		else
		{
			iter->pNode->cullChildren(cutoff);
			++iter;
		}
	}
}

// Return the number of nodes in the trie rooted at this node
size_t SeqTrie::Node::countNodes() const
{
	size_t count = 1;
	for(LinkList::const_iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
		count += iter->pNode->countNodes();
	return count;
}

// Recursive dot writer function
void SeqTrie::Node::writeDot(std::ostream& out) const
{
	out << "\"" << this << "\" [label=\"\"];\n";
	for(LinkList::const_iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
	{
		out << "\"" << this << "\" -> \"" << iter->pNode << "\" [label=\""
		    << iter->label << "," << iter->weight << "\"];\n";
		iter->pNode->writeDot(out);
	}
}

//
// SeqTrie
//

//
SeqTrie::SeqTrie()
{
	m_pRoot = new Node(NULL, '^');
}

//
SeqTrie::~SeqTrie()
{
	delete m_pRoot; // recursively destroys children
	m_pRoot = NULL;
}

// insert the string s into the trie
void SeqTrie::insert(const std::string& s, double weight)
{
	if(s.empty())
		return;
	m_pRoot->insert(s, weight, 0);
}

// remove s from the trie
void SeqTrie::remove(const std::string& s)
{
	WARN_ONCE("SeqTrie::remove only decrementing but not reaping");
	m_pRoot->remove(s, 0);
}

// return the number of nodes in the trie
size_t SeqTrie::countNodes() const
{
	return m_pRoot->countNodes();
}

// remove all sub-tries that have a link less than cutoff
void SeqTrie::cull(int cutoff)
{
	m_pRoot->cullChildren(cutoff);
}

// score string s against the trie
// place the result in out
void SeqTrie::score(const std::string& s, double p_error, PathScoreVector& out)
{
	PathScore start = {"", 0};
	m_pRoot->score(s, p_error, 0, start, out);
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
