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
#include <math.h>

//
// Link
//
SeqTrie::Link::Link() : pNode(NULL), label('\0'), count(0), weight(0.0f)
{

}

SeqTrie::Link::Link(Node* p, char l) : pNode(p), label(l), count(0), weight(0.0f)
{

}

void SeqTrie::Link::increment()
{
    ++count;
}

void SeqTrie::Link::decrement()
{
    --count;
}

void SeqTrie::Link::addWeight(double w)
{
    weight += w;
}

//
// Node
//
SeqTrie::Node::Node(Node* pParent, char parentLabel)
{
    parentLink.pNode = pParent;
    parentLink.label = parentLabel;
}

// Destructor, destroy the children of this node
SeqTrie::Node::~Node()
{
    for(LinkList::iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
        delete iter->pNode;
}

// Return a pointer to the link with label otherwise NULL
SeqTrie::Link* SeqTrie::Node::getLink(char label)
{
    for(LinkList::iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
        if(iter->label == label)
            return &(*iter);
    return NULL;
}

// Create a new child node
SeqTrie::Link* SeqTrie::Node::createChild(char label)
{
    Node* pChild = new Node(this, label);
    Link l(pChild, label);
    childLinks.push_back(l);
    return &childLinks.back();
}

// Score the string s against the trie, descending from this node
void SeqTrie::Node::score(const std::string& s, double lp_correct, 
                                                double lp_error, 
                                                double lp_missing,
                                                size_t idx, const PathScore& curr, PathScoreVector& out)
{
    if(s.size() == idx || childLinks.empty())
    {
        // Fill the rest of the path with missing node scores
        PathScore terminal = curr;
        terminal.path_corrected = terminal.path_sequence;
        for(size_t i = idx; i < s.size(); ++i)
        {
            terminal.path_corrected.append(1, s[i]);
            terminal.path_sequence.append(1, 'N');
            terminal.path_score += lp_missing;
            terminal.branch_score += lp_missing;
            terminal.probVector.push_back(lp_missing);
        }
        terminal.path_score += terminal.branch_score;
        out.push_back(terminal);
    }

    // Descend into each subtree
    for(LinkList::iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
    {
        PathScore branch = curr;
        if(iter->label == s[idx])
        {
            branch.path_score += lp_correct;
        }
        else
        {
            branch.path_score += lp_error;
        }

        branch.branch_score = log(1.0f - exp(iter->weight));

        branch.path_sequence.append(1, iter->label);
        ++branch.branch_length;
        branch.branch_cov += iter->count;
        branch.probVector.push_back(iter->weight);

        iter->pNode->score(s, lp_correct, lp_error, lp_missing, idx + 1, branch, out);
    }
}

// Get all the sequences into this subtrie and place them in svOut
void SeqTrie::Node::getSequences(std::string curr, StringVector& svOut) const
{
    // Reached a leaf, add the sequence to the vector
    if(childLinks.empty() && !curr.empty())
        svOut.push_back(curr);

    for(LinkList::const_iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
    {
        std::string child_str = curr + iter->label;
        iter->pNode->getSequences(child_str, svOut);
    }
}

// insert the string s into this node starting with the symbol at idx
// returns true if the sequence was successfully inserted
bool SeqTrie::Node::insert(const std::string& s, double weight, size_t idx)
{
    if(s.empty())
        return false;

    char b = s[idx];
    Link* pLink = getLink(b);
    if(pLink == NULL)
    {
        pLink = createChild(b);
    }
    
    pLink->increment();
    pLink->addWeight(weight);

    // Recurse
    if(++idx != s.size())
        return pLink->pNode->insert(s, weight, idx);
    else
        return true;
}

// remove the string s from the trie starting at pNode and character idx
bool SeqTrie::Node::remove(const std::string& s, size_t idx)
{
    char b = s[idx];
    Link* pLink = getLink(b);
    if(pLink != NULL)
    {
        pLink->decrement();
        if(++idx != s.size())
            return pLink->pNode->remove(s, idx);
    }
    return false;
}


// Remove children with count below cutoff
void SeqTrie::Node::cullChildren(int cutoff)
{
    LinkList::iterator iter = childLinks.begin(); 
    while(iter != childLinks.end())
    {
        if(childLinks.size() > 1 && iter->count < cutoff)
        {
            delete iter->pNode; // recursive
            iter = childLinks.erase(iter);
        }
        else
        {
            iter->pNode->cullChildren(cutoff);
            ++iter;
        }
    }
}

// Remodel the trie by remapping low-count children to high count branches
void SeqTrie::Node::remodel(int cutoff, double weight)
{
    // Split the link list into strong/weak links
    LinkList strongLinks;
    LinkList weakLinks;

    for(LinkList::iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
    {
        if(iter->count < cutoff)
            weakLinks.push_back(*iter);
        else
            strongLinks.push_back(*iter);
    }

    // For all the weak links, if there is a unique solid link, re-map the sequence of the weak link
    // trie into the solid branch, otherwise just remove it
    bool bRemodel = strongLinks.size() == 1;
    for(LinkList::iterator iter = weakLinks.begin(); iter != weakLinks.end(); ++iter)
    {
        if(bRemodel)
        {
            // Get the sequences in the weak trie
            StringVector sv;
            iter->pNode->getSequences("", sv);

            // Insert the strings into the strong trie
            for(size_t i = 0; i < sv.size(); ++i)
            {
                //std::cout << "Inserting: " << sv[i] << "\n";
                Link& strong = strongLinks.front();
                strong.pNode->insert(sv[i], weight, 0);
            }
        }

        // Destroy the weak trie
        delete iter->pNode;
    }

    childLinks = strongLinks;
    
    // Recurse
    for(LinkList::iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
        iter->pNode->remodel(cutoff, weight);
}

// Return the number of nodes in the trie rooted at this node
size_t SeqTrie::Node::countNodes() const
{
    size_t count = 1;
    for(LinkList::const_iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
        count += iter->pNode->countNodes();
    return count;
}

// Recursive dot writer function
void SeqTrie::Node::writeDot(std::ostream& out) const
{
    out << "\"" << this << "\" [label=\"\"];\n";
    for(LinkList::const_iterator iter = childLinks.begin(); iter != childLinks.end(); ++iter)
    {
        out << "\"" << this << "\" -> \"" << iter->pNode << "\" [label=\""
            << iter->label << "," << iter->weight << "\"];\n";
        iter->pNode->writeDot(out);
    }
}

//
// SeqTrie
//

//
SeqTrie::SeqTrie()
{
    m_pRoot = new Node(NULL, '^');
}

//
SeqTrie::~SeqTrie()
{
    delete m_pRoot; // recursively destroys children
    m_pRoot = NULL;
}

// insert the string s into the trie
void SeqTrie::insert(const std::string& s, double weight)
{
    if(s.empty())
        return;
    m_pRoot->insert(s, weight, 0);
}

// remove s from the trie
void SeqTrie::remove(const std::string& s)
{
    WARN_ONCE("SeqTrie::remove only decrementing but not reaping");
    m_pRoot->remove(s, 0);
}

// return the number of nodes in the trie
size_t SeqTrie::countNodes() const
{
    return m_pRoot->countNodes();
}

// remove all sub-tries that have a link less than cutoff
void SeqTrie::cull(int cutoff)
{
    m_pRoot->cullChildren(cutoff);
}

// remodel the trie but remapping low-probabilty branches
void SeqTrie::remodel(int cutoff, double weight)
{
    m_pRoot->remodel(cutoff, weight);
}

// score string s against the trie
// place the result in out
void SeqTrie::score(const std::string& s, double p_error, PathScoreVector& out)
{
    PathScore start;
    double lp_error = log(p_error);
    double lp_correct = log(1.0f - p_error);
    double lp_missing = log(0.9);
    m_pRoot->score(s, lp_correct, lp_error, lp_missing, 0, start, out);
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
