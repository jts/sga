//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//
// SearchHistory - data structures to track
// the history of a FM-index search
//
//-----------------------------------------------
#include "SearchHistory.h"
#include <algorithm>
#include <iterator>

// Link
SearchHistoryLink::SearchHistoryLink(SearchHistoryNode* ptr) : pNode(ptr)
{
	if(pNode != NULL)
	{
		pNode->increment();
	}
}

//
SearchHistoryLink::~SearchHistoryLink()
{
	if(pNode != NULL)
	{
		pNode->decrement();
		if(pNode->getCount() == 0)
			delete pNode;
	}
}

SearchHistoryLink::SearchHistoryLink(const SearchHistoryLink& link) : pNode(link.pNode)
{
	std::cout << "Copy CTOR\n";
	if(pNode != NULL)
		pNode->increment();
}

SearchHistoryLink& SearchHistoryLink::operator=(SearchHistoryLink& link)
{
	std::cout << "Assign op\n";
	SearchHistoryNode* pOld = pNode;
	pNode = link.pNode;
	if(pNode != NULL)
		pNode->increment();

	if(pOld != NULL)
	{
		pOld->decrement();
		if(pOld->getCount() == 0)
			delete pOld;
	}
	return *this;
}

//
// Node
//
SearchHistoryNode::SearchHistoryNode(SearchHistoryNode* pParent, 
                                     int var_pos, char var_base) : m_parentLink(pParent), 
									                               m_variant(var_pos, var_base),
																   m_refCount(0)
{
	
}


// 
SearchHistoryNode::~SearchHistoryNode()
{
	assert(m_refCount == 0);
	std::cout << "DELETE, node now: " << m_variant << " : " << m_refCount << "\n";
}

// adding a child of the node automatically increments the refCount of this node
// through the child's Link to this node. Once all the children of this node
// have been removed it will automatically be deleted
SearchHistoryLink SearchHistoryNode::createChild(int var_pos, char var_base)
{
	return SearchHistoryLink(new SearchHistoryNode(this, var_pos, var_base));
}

// Return the search history up to the root node
SearchHistoryVector SearchHistoryNode::getHistory()
{
	SearchHistoryVector out;
	SearchHistoryLink pCurr(this);
	bool done = false;
	while(!done)
	{
		if(pCurr->m_variant.base == '\0')
		{
			// root found, terminate the search
			done = true;
		}
		else
		{
			out.add(pCurr->m_variant);
			pCurr = pCurr->m_parentLink;
		}
	}
	return out;
}

//
void SearchHistoryNode::increment()
{
	++m_refCount;
	std::cout << "INC, node now: " << m_variant << " : " << m_refCount << "\n";

}

//
void SearchHistoryNode::decrement()
{
	assert(m_refCount != 0);
	--m_refCount;
	std::cout << "DEC, node now: " << m_variant << " : " << m_refCount << "\n";
}

// 
int SearchHistoryNode::getCount() const
{
	return m_refCount;
}

// SearchTree
SearchTree::SearchTree()
{
	m_pRoot = new SearchHistoryNode(NULL, -1, '\0');
	// The root should never be auto-destroyed so we set its count to one
	m_pRoot->increment();
}

//
SearchTree::~SearchTree()
{
	assert(m_pRoot->getCount() == 1);
	m_pRoot->decrement();
	delete m_pRoot;
}

//
SearchHistoryLink SearchTree::getRootLink()
{
	return SearchHistoryLink(m_pRoot);
}

//
// SearchHistoryVector
//

// Normalize the search history so that the positions are in ascending (left to right)
// order and the bases are wrt. the original strand of the forward sequence
// This is required to compare two SearchHistories that represent alignments to different 
// strands
void SearchHistoryVector::normalize(bool doComplement)
{
	std::sort(m_history.begin(), m_history.end(), SearchHistoryItem::sortPos);

	if(doComplement)
	{
		for(size_t i = 0; i < m_history.size(); ++i)
		{
			m_history[i].base = complement(m_history[i].base);	
		}
	}
}

//
void SearchHistoryVector::add(int pos, char base)
{
	SearchHistoryItem item(pos, base);
	m_history.push_back(item);
}

// 
void SearchHistoryVector::add(SearchHistoryItem& item)
{
	m_history.push_back(item);
}

//
int SearchHistoryVector::countDifferences(const SearchHistoryVector& a, const SearchHistoryVector& b, int maxPos)
{
	size_t na = a.m_history.size();
	size_t nb = b.m_history.size();

	size_t i = 0;
	size_t j = 0;
	int count = 0;

	while(i < na && j < nb)
	{
		const SearchHistoryItem& itemA = a.m_history[i];
		const SearchHistoryItem& itemB = b.m_history[j];
		
		if(itemA.pos > maxPos || itemB.pos > maxPos)
			break;

		if(itemA.pos == itemB.pos)
		{
			if(itemA.base != itemB.base)
				++count;
			++i;
			++j;
		}
		else if(itemA.pos < itemB.pos)
		{
			++count;
			++i;
		}
		else if(itemB.pos < itemA.pos)
		{
			++count;
			++j;
		}
		else
		{
			assert(false);
		}
	}
	
	// Count the remaining elements in A and B that are less than maxpos
	while(i < na && a.m_history[i].pos <= maxPos)
	{
		++count;
		++i;
	}

	// Count the remaining elements in A and B that are less than maxpos
	while(j < nb && b.m_history[j].pos <= maxPos)
	{
		++count;
		++j;
	}

	return count;
}

//
std::ostream& operator<<(std::ostream& out, const SearchHistoryVector& hist)
{
	std::copy(hist.m_history.begin(), hist.m_history.end(), 
	          std::ostream_iterator<SearchHistoryItem>(out, "\t"));
	return out;
}

