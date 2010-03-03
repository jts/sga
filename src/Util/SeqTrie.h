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

struct PathScore
{
	std::string path;
	double score;
};
typedef std::vector<PathScore> PathScoreVector;


class SeqTrie
{
	// Internal datastructures
	struct Node;
	struct Link
	{
		// functions
		Link();
		Link(Node* p, char l);
		void increment();
		void decrement();
		void addWeight(double w);

		// data
		Node* pNode;
		char label;
		int count;
		double weight;
	};

	typedef std::list<Link> LinkList;

	class Node
	{
		public:
			// functions
			Node(Node* pParent, char parentLabel);
			~Node();

			Link* getLink(char label);

			bool insert(const std::string& s, double weight, size_t idx);
			bool insertAtDepth(const std::string& s, double weight, size_t depth);
			bool remove(const std::string& s, size_t idx);
			
			void score(const std::string& s, double p_error, size_t idx, const PathScore& curr, PathScoreVector& out);

			void cullChildren(int cutoff);
			void writeDot(std::ostream& out) const;

		private:
	
			Link* createChild(char label);
				
			//data
			Link parentLink;
			LinkList pChildLinks;
	};
	
	// 
	public:

		SeqTrie();
		~SeqTrie();

		void score(const std::string& s, double p_error, PathScoreVector& out);
	
		// Creation functions
		void insert(const std::string& s, double weight);
		void insertAtDepth(const std::string& s, double weight, size_t depth);
	
		// Remove the string s from the trie
		void remove(const std::string& s);
		void cull(int cutoff);

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

