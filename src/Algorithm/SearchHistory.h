//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//
// SearchHistory - data structures to track
// the history of a FM-index search
//-----------------------------------------------
#ifndef SEARCHHISTORY_H
#define SEARCHHISTORY_H

#include "Util.h"

//
struct SearchHistoryItem
{
	SearchHistoryItem(int p, char b) : pos(p), base(b) {}
	int pos;
	char base;

	static bool sortPos(const SearchHistoryItem& a, const SearchHistoryItem& b)
	{
		return a.pos < b.pos;
	}

	friend std::ostream& operator<<(std::ostream& out, const SearchHistoryItem& hi)
	{
		out << hi.pos << "," << hi.base;
		return out;
	}
};
typedef std::vector<SearchHistoryItem> HistoryItemVector;

//
class SearchHistoryVector
{
	public:
		
		//
		SearchHistoryVector() {}

		//
		void add(int pos, char base);
		void add(SearchHistoryItem& item);
		void normalize(bool doComplement);

		//
		static int countDifferences(const SearchHistoryVector& a, const SearchHistoryVector& b, int maxPos);

		//
		friend std::ostream& operator<<(std::ostream& out, const SearchHistoryVector& hist);

	private:

		HistoryItemVector m_history;
};


class SearchHistoryNode;
typedef std::vector<SearchHistoryNode*> NodeVector;

// A SearchHistoryLink is a reference-counted wrapper of a 
// search node. This is used as the interface to the search tree
// so that nodes are deleted when they are no longer used.
class SearchHistoryLink
{
	public:
		SearchHistoryLink(SearchHistoryNode* ptr);

		// To get the reference counting correct we need to handle the copy 
		// constructor and the assignment operator
		SearchHistoryLink(const SearchHistoryLink& link);
		SearchHistoryLink& operator=(SearchHistoryLink& link);
		~SearchHistoryLink();

		SearchHistoryNode* operator-> () { assert(pNode != NULL); return pNode; }
		SearchHistoryNode& operator* ()  { assert(pNode != NULL); return *pNode; }

	private:
		SearchHistoryLink() {} // Not allowed
		SearchHistoryNode* pNode;
};

//
class SearchHistoryNode
{
	public:
		SearchHistoryLink createChild(int var_pos, char var_base);
		SearchHistoryVector getHistory();

	private:

		friend class SearchHistoryLink;
		friend class SearchTree;

		// The nodes should only be constructed/destructed through the SearchTree Links
		SearchHistoryNode(SearchHistoryNode* pParent, int var_pos, char var_base);
		~SearchHistoryNode();

		void increment();
		void decrement();
		int getCount() const;

		SearchHistoryLink m_parentLink;
		SearchHistoryItem m_variant;
		int m_refCount;
};

//
class SearchTree
{
	public:

		SearchTree();
		~SearchTree();

		SearchHistoryLink getRootLink();

	private:
		SearchHistoryNode* m_pRoot;
};


#endif
