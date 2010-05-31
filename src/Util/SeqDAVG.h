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
        class Node;
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
                Link* addLink(Node* pNode, double weight, char label);
                void writeDot(std::ostream& out) const;

            private:
                    
                //data
                LinkList pChildLinks;
        };
    
        typedef std::list<Node*> NodePList;

        //
        SeqDAVG();
        ~SeqDAVG();

        //
        void insert(const std::string& s, double weight, size_t depth = 0);
        
        // I/O
        void writeDot(std::string filename);


    private:
        static Link* find(LinkList& list, char label);
        static LinkList findList(LinkList& list, char label);

        std::vector<LinkList> m_data;
        Node* m_pRoot;
};

#endif
