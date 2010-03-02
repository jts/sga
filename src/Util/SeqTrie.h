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
#include <list>
class SeqTrie
{
	// Internal datastructures
	struct Node;
	struct Link
	{
		// functions
		void increment();
		void decrement();

		// data
		Node* pNode;
		char label;
		int count;
	};
	
	typedef std::list<Link> LinkList;

	struct Node
	{
		// functions
		Node(Node* pParent, char parentLabel);
		~Node();

		Link* getLink(char label);
		Link* createChild(char label);
		void removeChildren(int cutoff);
		void writeDot(std::ostream& out) const;

		//data
		Link parentLink;
		LinkList pChildLinks;
	};
	
	// 
	public:

		SeqTrie();
		~SeqTrie();

		// Creation functions
		void insert(const std::string& s);
		void insertAtDepth(const std::string& s, size_t depth);
		
		// Remove the string s from the trie
		void remove(const std::string& s);
		void reap(int cutoff);

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

