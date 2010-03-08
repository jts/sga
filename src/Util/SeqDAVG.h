//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SeqDAVG.h - Directed acyclic variant graph
// for a fixed-length sequence. Contains all the variants
// seen in the overlaps for a single string.
//
#ifndef SEQDAVG_H
#define SEQDAVG_H

#include "Util.h"
#include <list>

class SeqDAVG
{
	public:
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
				Node();
				~Node();

				Link* getLink(char label);

				bool insert(const std::string& s, double weight, size_t idx);
				//bool remove(const std::string& s, size_t idx);
				//void score(const std::string& s, double p_error, size_t idx, const PathScore& curr, PathScoreVector& out);

				void cullChildren(int cutoff);
				void writeDot(std::ostream& out) const;

			private:
		
				Link* createChild(char label);
					
				//data
				LinkList pChildLinks;
		};


		//
		SeqDAVG(const int size);
		~SeqDAVG();

		//
		void insert(const std::string& s);
		void insertAtDepth(const std::string& s);

		void cull(int cutoff);
		
		// I/O
		void writeDot(std::string filename);


	private:

		std::vector<LinkList> m_data;

};

#endif
