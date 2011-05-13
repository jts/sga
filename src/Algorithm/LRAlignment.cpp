//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// LRAlignment - Collection of algorithms for performing
// long read matches against an FM-index
//
#include "LRAlignment.h"
#include "QuickBWT.h"
namespace LRAlignment
{
static const int MINUS_INF = -0x3fffffff;

// Initialize the LRCell with default values
void LRCell::initializeDefault()
{
    interval.lower = 0;
    interval.upper = 0;
    I = MINUS_INF;
    D = MINUS_INF;
    G = MINUS_INF;

    parent_cidx = 0;
    q_len = 0;
    t_len = 0;
    parent_idx = -1;
    u_idx = -1;

    clearChildren();
}

void LRCell::clearChildren()
{
    children_idx[0] = -1;
    children_idx[1] = -1;
    children_idx[2] = -1;
    children_idx[3] = -1;
}

bool LRCell::hasUninitializedChild() const
{
    return children_idx[0] == -1 || children_idx[1] == -1 || children_idx[2] == -1 ||  children_idx[3] == -1;
}

// Implementation of bwa-sw algorithm.
// Roughly follows Heng Li's code
void bwaswAlignment(const std::string& query, BWT* pTargetBWT)
{
    // Construct an FM-index of the query sequence
    BWT* pQueryBWT = createQuickBWT(query);

    LRParams params;

    // Initialize a stack of elements with a single entry for the root node
    // of the query DAWG
    LRStack stack;
    
    //
    LRStackEntry initial;
    initial.interval.lower = 0;
    initial.interval.upper = pQueryBWT->getBWLen() - 1;

    //
    LRCell x;
    x.initializeDefault();
    x.G = 0;
    x.interval.lower = 0;
    x.interval.upper = pTargetBWT->getBWLen() - 1;
    x.q_len = pTargetBWT->getBWLen();
    initial.cells.push_back(x);

    stack.push(initial);

    int abort = 5;
    // Traverse DAG
    while(!stack.empty())
    {
        LRStackEntry v = stack.top();
        stack.pop();
        
        size_t old_n = v.cells.size();
        printf("Outer stack loop [%d %d] old_n: %d\n", old_n, v.interval.lower, v.interval.upper);

        // TODO: bandwidth test and max depth
        if(abort-- == 0)
            exit(1);
        
        // Calculate occurrence values for children of the current entry
        for(int qci = 0; qci < DNA_ALPHABET::size; ++qci)
        {
            char query_child_base = DNA_ALPHABET::getBase(qci);
            BWTInterval child_interval = v.interval;
            BWTAlgorithms::updateInterval(child_interval, query_child_base, pQueryBWT);
            std::cout << "Query child interval(" << query_child_base << "): " << child_interval << "\n";
            if(!child_interval.isValid())
                continue;

            // TODO: hash stuff?

            // Create an array of cells for the current node holding
            // scores between the current node and nodes in the prefix trie
            LRStackEntry u;
            u.interval = child_interval;

            // Loop over the nodes in v
            for(size_t i = 0; i < v.cells.size(); ++i)
            {
                LRCell* p = &v.cells[i];
                LRCell x; // the cell being calculated
                LRCell* c[4]; // pointers to cells required to calculate x

                c[0] = &x;
                x.G = MINUS_INF;
                x.u_idx = p->u_idx = -1;

                bool add_x_to_u = false;
                if(p->parent_idx >= 0) 
                {
                    std::cout << "parent has been visited\n";

                    int upos = v.cells[p->parent_idx].u_idx;
                    std::cout << "  upos: " << upos << "\n";
                    c[1] = upos >= 0 ? &u.cells[upos] : NULL;
                    c[2] = p;
                    c[3] = &v.cells[p->parent_idx];

                    int match_score = qci == p->parent_cidx ? params.match : -params.mismatch;
                    std::cout << "  qci: " << qci << " pci: " << (int)p->parent_cidx << "\n";
                    std::cout << "  match score: " << match_score << "\n";
                    int score = fill_cells(params, match_score, c);
                    std::cout << "  score: " << score << "\n";
                    if(score > 0)
                    {
                        // this cell has a positive score
                        x.parent_idx = upos;
                        p->u_idx = u.cells.size(); // x will be added to u in this position
                        add_x_to_u = true;
                    }
                }
                else
                {
                    std::cout << "parent does not exist\n";
                    if(p->D > p->G - params.gap_open)
                        x.D = p->D - params.gap_ext; // extend gap
                    else
                        x.D = p->G - params.gap_open_extend; // open a new gap

                    if(x.D > 0)
                    {
                        x.G = x.D;
                        x.I = MINUS_INF;
                        x.parent_idx = -1;
                        p->u_idx = u.cells.size();
                        add_x_to_u = true;
                    }
                }

                if(add_x_to_u)
                {
                    std::cout << "Adding cell\n";
                    // set the remaining fields in the current cell
                    x.clearChildren();
                    x.parent_cidx = p->parent_cidx;
                    x.interval = p->interval;
                    x.q_len = p->q_len + 1;
                    x.t_len = p->t_len;
                    u.cells.push_back(x);

                    // TODO: some heap stuff?
                }

                // Check if we should descend into another node of the prefix trie of the target
                std::cout << "X.G: " << x.G << " i: " << i << " oldn: " << old_n << "\n";
                if( (x.G > params.gap_open_extend /* && heap test */) || i < old_n)
                {
                    if(p->hasUninitializedChild())
                    {
                        for(int tci = 0; tci < DNA_ALPHABET::size; ++tci)
                        {
                            if(p->children_idx[tci] != -1)
                                continue; // already added
                            char target_child_base = DNA_ALPHABET::getBase(tci);
                            BWTInterval target_child_interval = p->interval;
                            BWTAlgorithms::updateInterval(target_child_interval, target_child_base, pTargetBWT);
                            std::cout << "Target child interval(" << target_child_base << "): " << target_child_interval << "\n";
                            if(!target_child_interval.isValid()) // child with this extension does not exist
                            {
                                p->children_idx[tci] = -2;
                                continue;
                            }

                            // Create new entry on array v
                            LRCell y;
                            y.G = y.I = y.D = MINUS_INF;
                            y.interval = target_child_interval;
                            y.parent_cidx = tci;
                            y.q_len = p->q_len;
                            y.t_len = p->t_len + 1;
                            y.parent_idx = i;
                            y.clearChildren();
                            p->children_idx[tci] = v.cells.size();

                            v.cells.push_back(y);

                            p = &v.cells[i]; // p may have been invalidated in the push, update
                        } // for tcl
                    } // if has uninitialized
                } // if X.h
            } // for all v
            // push u onto the stack
            stack.push(u);
            printf("pushing u to some stack [%d %d]\n", u.interval.lower, u.interval.upper);


        } // for qci
    } // for all stack
    delete pQueryBWT;
}

// Fill the values of C[0] depending on the values in the other 3 cells
int fill_cells(const LRParams& params, int match_score, LRCell* c[4])
{
	int G = c[3] ? c[3]->G + match_score : MINUS_INF;
	if(c[1]) 
    {
		if(c[1]->I > c[1]->G - params.gap_open)
            c[0]->I = c[1]->I - params.gap_ext; // extend gap
        else
            c[0]->I = c[1]->G - params.gap_open_extend; // open new gap
		if (c[0]->I > G)
            G = c[0]->I; // new best score
	} 
    else
    {
        c[0]->I = MINUS_INF;
    }

	if(c[2])
    {
		if(c[2]->D > c[2]->G - params.gap_open)
            c[0]->D = c[2]->D - params.gap_ext; // extend gap
        else
            c[0]->D = c[2]->G - params.gap_open_extend; // open new gap

		if (c[0]->D > G) 
            G = c[0]->D; // new best score
	} 
    else
    {
        c[0]->D = MINUS_INF;
    }
    
	return(c[0]->G = G);
}

};
