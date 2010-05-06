//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//
// SearchTree - data structure to track
// the history of a FM-index search
//
//-----------------------------------------------
#include "SearchTree.h"

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

// Node
SearchHistoryNode::SearchHistoryNode(SearchHistoryNode* pParent, 
                                     int var_pos, char var_base) : m_parentLink(pParent), 
									                               m_variant(var_pos, var_base),
																   m_refCount(0)
{
	
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
SearchHistoryNode::~SearchHistoryNode()
{
	assert(m_refCount == 0);
	std::cout << "DELETE, node now: " << m_variant << " : " << m_refCount << "\n";
}

// adding a child of the node automatically increments the refCount of this node
// through the child's Link to this node. Once all the children of this node
// have been removed it will automatically be deleted
SearchHistoryLink SearchHistoryNode::addChild(int var_pos, char var_base)
{
	return SearchHistoryLink(new SearchHistoryNode(this, var_pos, var_base));
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
	m_pRoot = new SearchHistoryNode(NULL, -1, 'X');
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
