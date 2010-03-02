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

class SeqTrie
{
	// Internal datastructures
	struct Node;
	struct Link
	{
		// functions
		void addWeight(double w);

		// data
		Node* pNode;
		char label;
		double weight;
	};

	typedef std::vector<Link> LinkVector;

	struct Node
	{
		// functions
		Node(Node* pParent, char parentLabel);
		~Node();

		Link* getLink(char label);
		Link* createChild(char label);
		void writeDot(std::ostream& out) const;

		//data
		Link parentLink;
		std::vector<Link> pChildLinks;
	};
	
	// 
	public:

		SeqTrie();
		~SeqTrie();

		// Creation functions
		void insert(const std::string& s);

		// I/O
		void writeDot(std::string filename);

	private:

		void insert(Node* pNode, const std::string& s, size_t idx);

		Node* m_pRoot;
};

#endif

