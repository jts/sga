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
void SeqTrie::Link::addWeight(double w)
{
	weight += w;
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
	for(size_t i = 0; i < pChildLinks.size(); ++i)
		delete pChildLinks[i].pNode;
}

// Return a pointer to the link with label otherwise NULL
SeqTrie::Link* SeqTrie::Node::getLink(char label)
{
	for(size_t i = 0; i < pChildLinks.size(); ++i)
	{
		if(pChildLinks[i].label == label)
			return &pChildLinks[i];
	}
	return NULL;
}

// Create a new child node
SeqTrie::Link* SeqTrie::Node::createChild(char label)
{
	Node* pChild = new Node(this, label);
	Link l = {pChild, label, 0.0f};
	pChildLinks.push_back(l);
	return &pChildLinks.back();
}

// Recursive dot writer function
void SeqTrie::Node::writeDot(std::ostream& out) const
{
	out << "\"" << this << "\" [label=\"\"];\n";
	for(size_t i = 0; i < pChildLinks.size(); ++i)
	{
		out << "\"" << this << "\" -> \"" << pChildLinks[i].pNode << "\" [label=\""
		    << pChildLinks[i].label << "," << pChildLinks[i].weight << "\"];\n";
		pChildLinks[i].pNode->writeDot(out);
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

// insert the string s into the trie starting at pNode and character idx
void SeqTrie::insert(Node* pNode, const std::string& s, size_t idx)
{
	assert(pNode != NULL);
	char b = s[idx];
	Link* pLink = pNode->getLink(b);
	if(pLink == NULL)
	{
		// The node does not have a link with this label, create a new node
		pLink = pNode->createChild(b);
	}
	pLink->addWeight(1.0f);
	// Recurse
	if(++idx != s.size())
		return insert(pLink->pNode, s, idx);
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
