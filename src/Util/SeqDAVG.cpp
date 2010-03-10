//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
//
// SeqDAVG.h - Directed acyclic variant graph
// for a fixed-length sequence. Contains all the variants
// seen in the overlaps for a single string.
//
#include "SeqDAVG.h"
#include <iostream>

//
// Link
//
SeqDAVG::Link::Link() : pNode(NULL), label('\0'), count(0), weight(0.0f)
{

}

SeqDAVG::Link::Link(Node* p, char l) : pNode(p), label(l), count(0), weight(0.0f)
{

}

void SeqDAVG::Link::increment()
{
	++count;
}

void SeqDAVG::Link::decrement()
{
	--count;
}

void SeqDAVG::Link::addWeight(double w)
{
	weight += w;
}

//
// Node
//
SeqDAVG::Node::Node()
{
}

//
SeqDAVG::Node::~Node()
{
}

// Return a pointer to the link with label otherwise NULL
SeqDAVG::Link* SeqDAVG::Node::getLink(char label)
{
	return SeqDAVG::find(pChildLinks, label);
}

// Create a new child node
SeqDAVG::Link* SeqDAVG::Node::addLink(Node* pNode, double weight, char label)
{
	Link* pLink = getLink(label);
	if(pLink == NULL)
	{
		pChildLinks.push_back(Link(pNode, label));
		pLink = &pChildLinks.back();
	}
	pLink->addWeight(weight);
	pLink->increment();
	return pLink;
}

// Output node in dot format
void SeqDAVG::Node::writeDot(std::ostream& out) const
{
	out << "\"" << this << "\" [label=\"\"];\n";
	for(LinkList::const_iterator iter = pChildLinks.begin(); iter != pChildLinks.end(); ++iter)
	{
		out << "\"" << this << "\" -> \"" << iter->pNode << "\" [label=\""
		    << iter->label << "," << iter->weight << "\"];\n";
	}
}

//
// SeqDAVG
//

//
SeqDAVG::SeqDAVG()
{
	m_pRoot = new Node();
}

//
SeqDAVG::~SeqDAVG()
{
	delete m_pRoot;
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		for(LinkList::iterator iter = m_data[i].begin(); iter != m_data[i].end(); ++iter)
			delete iter->pNode;
	}
}

// insert the string s into the trie so that it is a child 
// of the node(s) at DEPTH. Children of the root (the first
// nodes) are depth 0. If depth is higher than the deepest
// node in the trie, this will do nothing.
void SeqDAVG::insert(const std::string& s, double weight, size_t depth)
{
	if(s.empty())
		return;

	// Expand the data vector if necessary
	if(depth + s.size() > m_data.size())
	{
		m_data.resize(depth + s.size());
	}

	char prev_label = '*';

	for(size_t i = 0; i < s.size(); ++i)
	{
		size_t curr_depth = depth + i;
		char label = s[i];

		// Find the node, if it doesnt exist create it
		Link* pNodeLink = find(m_data[curr_depth], label);
		if(pNodeLink == NULL)
		{
			// Create the new node
			Node* pNode = new Node;
			m_data[curr_depth].push_back(Link(pNode, label));
			pNodeLink = &m_data[curr_depth].back();
		}

		// Link this node to previous nodes - the '*' label matches everything. 
		// If depth is zero just link to the root
		if(curr_depth == 0)
		{
			m_pRoot->addLink(pNodeLink->pNode, weight, pNodeLink->label);
		}
		else
		{
			LinkList& list = m_data[curr_depth - 1];
			for(LinkList::iterator iter = list.begin(); iter != list.end(); ++iter)
				if(prev_label == '*' || iter->label == prev_label)
					iter->pNode->addLink(pNodeLink->pNode, weight, pNodeLink->label);
		}
		prev_label = label;
	}
}

//
SeqDAVG::Link* SeqDAVG::find(LinkList& list, char label)
{
	for(LinkList::iterator iter = list.begin(); iter != list.end(); ++iter)
		if(iter->label == label)
			return &(*iter);
	return NULL;

}

// Write the trie to a dot file
void SeqDAVG::writeDot(std::string filename)
{
	std::ofstream writer(filename.c_str());
	writer << "digraph G\n{\n";
	m_pRoot->writeDot(writer);
	for(size_t i = 0; i < m_data.size(); ++i)
	{
		for(LinkList::iterator iter = m_data[i].begin(); iter != m_data[i].end(); ++iter)
			iter->pNode->writeDot(writer);
	}
	writer << "}\n";
	writer.close();
}
