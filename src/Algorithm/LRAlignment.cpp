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
    BWT* pQueryBWT = NULL;
    SuffixArray* pQuerySA = NULL;
    createQuickBWT(query, pQueryBWT, pQuerySA);

    LRHash dupHash;

    // Initialize the hash table of DAWG nodes
    LRHash dawgHash;
    initializeDAWGHash(pQueryBWT, dawgHash);
    
    //
    LRParams params;

    // Initialize a stack of elements with a single entry for the root node
    // of the query DAWG
    LRStack stack;
    
    // Initialize hits vector with enough space to hold 2 hits per query base
    LRHitVector hitsVector(2*pQueryBWT->getBWLen());

    // Initialize the vector of dawg nodes pending addition to the stack
    LRPendingVector pendingVector;
    size_t num_pending = 0;

    //
    LRStackEntry* pInitial = new LRStackEntry;
    pInitial->interval.lower = 0;
    pInitial->interval.upper = pQueryBWT->getBWLen() - 1;

    //
    LRCell x;
    x.initializeDefault();
    x.G = 0;
    x.interval.lower = 0;
    x.interval.upper = pTargetBWT->getBWLen() - 1;
    x.t_len = 0;
    x.q_len = 0;

    pInitial->cells.push_back(x);

    stack.push(pInitial);

    // Traverse DAG
    while(!stack.empty())
    {
        LRStackEntry* v = stack.top();
        stack.pop();
        
        size_t old_n = v->cells.size();
        printf("Outer stack loop [%zu %zu] old_n: %zu\n", old_n, v->interval.lower, v->interval.upper);

        // TODO: bandwidth test and max depth
        
        // Calculate occurrence values for children of the current entry
        for(int qci = 0; qci < DNA_ALPHABET::size; ++qci)
        {
            char query_child_base = DNA_ALPHABET::getBase(qci);
            BWTInterval child_interval = v->interval;
            BWTAlgorithms::updateInterval(child_interval, query_child_base, pQueryBWT);
            std::cout << "Query child interval(" << qci << "): [" << child_interval << "]\n";
            if(!child_interval.isValid())
                continue;

            // Update the hashed count of the number of times this substring needs to be visited
            uint64_t key = child_interval.lower << 32 | child_interval.upper;
            LRHash::iterator hashIter = dawgHash.find(key);
            assert(hashIter != dawgHash.end());
            assert((uint32_t)hashIter->second > 0);
            --hashIter->second;

			printf("hash k1: %zu k2: %zu v1: %zu v2: %d\n", child_interval.lower,  child_interval.upper, hashIter->second >> 32, (uint32_t)hashIter->second);

            // Create an array of cells for the current node holding
            // scores between the current node and nodes in the prefix trie
            LRStackEntry* u = new LRStackEntry;
            u->interval = child_interval;

            // Loop over the nodes in v
            for(size_t i = 0; i < v->cells.size(); ++i)
            {
                LRCell* p = &v->cells[i];
                printf("processing p = v[%zu] [%zu %zu] is deleted? %d\n", i, p->interval.lower == -1 ? 0 : p->interval.lower, p->interval.upper == -1 ? 0 : p->interval.upper, p->interval.upper == -1);
                if(p->interval.upper == -1)
                    continue; // duplicate that has been deleted
                LRCell x; // the cell being calculated
                LRCell* c[4]; // pointers to cells required to calculate x

                c[0] = &x;
                x.G = MINUS_INF;
                x.u_idx = p->u_idx = -1;

                bool add_x_to_u = false;
                if(p->parent_idx >= 0) 
                {
                    std::cout << "parent has been visited\n";

                    int upos = v->cells[p->parent_idx].u_idx;
                    std::cout << "    upos: " << upos << "\n";
                    c[1] = upos >= 0 ? &u->cells[upos] : NULL;
                    c[2] = p;
                    c[3] = &v->cells[p->parent_idx];

                    int match_score = qci == p->parent_cidx ? params.match : -params.mismatch;
                    int score = fillCells(params, match_score, c);
                    std::cout << "    score: " << score << "\n";
                    if(score > 0)
                    {
                        // this cell has a positive score
                        // set the x's parent position in u and set the
                        // parent's child position to x
                        x.parent_idx = upos;
                        p->u_idx = u->cells.size(); // x will be added to u in this position
                        
                        if(x.parent_idx >= 0)
                            u->cells[upos].children_idx[p->parent_cidx] = p->u_idx;
                        add_x_to_u = true;
                    }
                }
                else
                {
                    std::cout << "parent not visited\n";
                    if(p->D > p->G - params.gap_open)
                        x.D = p->D - params.gap_ext; // extend gap
                    else
                        x.D = p->G - params.gap_open_extend; // open a new gap

                    if(x.D > 0)
                    {
                        x.G = x.D;
                        x.I = MINUS_INF;
                        x.parent_idx = -1;
                        p->u_idx = u->cells.size();
                        add_x_to_u = true;
                    }
                }

                if(add_x_to_u)
                {
                    std::cout << "adding x to u\n";
                    // set the remaining fields in the current cell
                    x.clearChildren();
                    x.parent_cidx = p->parent_cidx;
                    x.interval = p->interval;
                    x.q_len = p->q_len + 1;
                    x.t_len = p->t_len;
                    u->cells.push_back(x);

                    // TODO: some heap stuff?
                }

                // Check if we should descend into another node of the prefix trie of the target
                std::cout << "x->G: " << x.G << " qr: " << params.gap_open_extend << " i: " << i << " old_n: " << old_n << "\n";
                if( (x.G > params.gap_open_extend /* && heap test */) || i < old_n)
                {
                    printf("    pass test 1\n");
                    if(p->hasUninitializedChild())
                    {
                        printf("    pass test 2\n");
                        for(int tci = 0; tci < DNA_ALPHABET::size; ++tci)
                        {
                            if(p->children_idx[tci] != -1)
                            {
                                printf("    child exists at pos %d\n", p->children_idx[tci]);
                                continue; // already added
                            }
                            else
                            {
                                printf("    child does not exist %d\n", p->children_idx[tci]);
                            }
                            char target_child_base = DNA_ALPHABET::getBase(tci);
                            BWTInterval target_child_interval = p->interval;
                            BWTAlgorithms::updateInterval(target_child_interval, target_child_base, pTargetBWT);
                            std::cout << "Target child interval(" << tci << "): [" << target_child_interval << "]\n";
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

                            p->children_idx[tci] = v->cells.size();
                            v->cells.push_back(y);

                            p = &v->cells[i]; // p may have been invalidated in the push, update
                        } // for tcl
                    } // if has uninitialized
                } // if X.h
            } // for all v
		
            // If the score for any cell in u exceeds the threshold, save it
            if(!u->cells.empty())
                saveHits(pQuerySA, u, params.threshold, hitsVector);

            //
            uint32_t count = (uint32_t)hashIter->second;
            uint32_t position = hashIter->second >> 32;
            printf("cnt: %d pos: %d\n", count, position);

            // Check if an entry in the pending array exists for this DAWG node
            if(position > 0)
            {
                LRStackEntry* w = pendingVector[position - 1];

                // An entry in the pending array has been created for this query substing
                // Merge u into the interval
                if(!u->cells.empty())
                {
                    // if the pending value has fewer cells than u
                    // swap their pointers
                    if(w->cells.size() < u->cells.size())
                    {
                        w = u;
                        u = pendingVector[position - 1];
                        pendingVector[position - 1] = w;
                    }

                    printf("merging stack entries\n");
                    mergeStackEntries(w, u);
                }

                if(count == 0)
                {
                    // this node in the dawg will not be visited again
                    // move the stack entry from the pending list to the stack
                    removeDuplicateCells(w, dupHash);
                    printf("moving w from pending to stack[%zu %zu]\n", w->interval.lower, w->interval.upper);
                    stack.push(w);
                    pendingVector[position - 1] = 0;
                    num_pending -= 1;
                }

                // u is empty or merged, it is no longer needed
                delete u;

            }
            else if(count > 0)
            {
                // Create an entry in the pending queue for the current node of the DAWG
                if(!u->cells.empty())
                {
                    pendingVector.push_back(u);
                    num_pending += 1;

                    // Save the position of u in the pending vector into the hash
                    // all subsequent traversals that visit this node of the dawg
                    // will get merged into this position. index + 1 is stored
                    // so that position == 0 indicates the empty case
                    hashIter->second = (uint64_t)pendingVector.size() << 32 | count;
                    printf("saving u to pending [%zu %zu]\n", u->interval.lower, u->interval.upper);
                }
                else
                {
                    // u has no cells to calculate, discard it
                    delete u;
                }
            }
            else // count == 0, pos == 0
            {
                // This substring is unique, push u straight onto the stack
                stack.push(u);
                printf("pushing u to stack [%zu %zu]\n", u->interval.lower, u->interval.upper);
            }
        } // for qci
        
        // done with v
        delete v;
    } // for all stack

    assert(num_pending == 0);

    delete pQueryBWT;
    delete pQuerySA;

    WARN_ONCE("CHECK FOR MEM LEAKS");
}

// Process a list of cells and save alignment hits meeting the threshold
void saveHits(const SuffixArray* pQuerySA, LRStackEntry* u, int threshold, LRHitVector& hits)
{
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        LRCell* p = &u->cells[i];
        if(p->G < threshold)
            continue;
        
        for(int64_t k = u->interval.lower; k <= u->interval.upper; ++k)
        {
            // Calculate the beginning of the alignment using the suffix array
            // of the query
            SAElem e = pQuerySA->get(k);
            int beg = e.getPos();
            int end = beg + p->q_len;

            // Save the best hit for alignments starting at beg in positions hits[2*beg]
            // and the second best hit in hits[2*beg+1]
            LRHit* q = NULL;
            if(p->G > hits[beg*2].G)
            {
                // move the previous best to the second best slot
                hits[beg*2+1] = hits[beg*2];
                q = &hits[2*beg];
            }
            else if(p->G > hits[beg*2+1].G)
            {
                // the current hit is the second best at this position
                q = &hits[2*beg+1];
            }

            if(q)
            {
                q->interval = p->interval;
                q->length = p->t_len;
                q->G = p->G;
                q->beg = beg;
                q->end = end;
                q->G2 = (q->interval.size() == 1) ? 0 : q->G;
                q->flag = 0;
                q->num_seeds = 0;
                printf("saved hit (%d %d) len: %d score: %d\n", q->beg, q->end, q->length, q->G);
            }
        }
    }
}

// Initialize the hash table of DAWG nodes by inserting an element
// into the hash table for each distinct substring of the query
void initializeDAWGHash(BWT* pQueryBWT, LRHash& hashTable)
{
    std::stack<uint64_t> stack;
    stack.push(pQueryBWT->getBWLen() - 1);
    while(!stack.empty())
    {
        uint64_t x = stack.top();
        stack.pop();

        BWTInterval interval(x >> 32, (uint32_t)x);
        AlphaCount64 lower = pQueryBWT->getFullOcc(interval.lower - 1);
        AlphaCount64 upper = pQueryBWT->getFullOcc(interval.upper);
        for(int ci = 0; ci < DNA_ALPHABET::size; ++ci)
        {
            char b = DNA_ALPHABET::getBase(ci);    
            size_t pb = pQueryBWT->getPC(b);
            interval.lower = pb + lower.get(b);
            interval.upper = pb + upper.get(b) - 1;

            if(!interval.isValid())
                continue;
            x = (uint64_t)interval.lower << 32 | interval.upper;

            // insert (or update) this key in the hash
            LRHash::iterator iter = hashTable.find(x);
            if(iter == hashTable.end())
            {
                hashTable.insert(std::make_pair(x,1));
                stack.push(x);
            }
            else
            {
                iter->second++;
            }
        }
    }
}

// Merge two stack entries which represent the same node in the DAWG
void mergeStackEntries(LRStackEntry* u, LRStackEntry* v)
{
    // Update parent index and child entries for v cell's to reflect
    // their new positions at the end of array of u
    size_t num_u = u->cells.size();
	for (size_t i = 0; i != v->cells.size(); ++i) 
    {
        LRCell* p = &v->cells[i];
		if(p->parent_idx >= 0)
            p->parent_idx += num_u;
        for(size_t j = 0; j < 4; ++j)
    		if (p->children_idx[j] >= 0) 
                p->children_idx[j] += num_u;
	}
    u->cells.insert(u->cells.end(), v->cells.begin(), v->cells.end());
}

// Remove duplicated cells from the stack entry
void removeDuplicateCells(LRStackEntry* u, LRHash& hash)
{
    hash.clear();
    int n = 0;
    for(size_t i = 0; i < u->cells.size(); ++i)
    {
        int j = -1; // element to be deleted
        LRCell* p = &u->cells[i];

        if(p->interval.upper == -1)
            continue; // already deleted

        uint64_t key = p->interval.lower << 32 | p->interval.upper;
        LRHash::iterator iter = hash.find(key);
        if(iter != hash.end())
        {
            // Check if the score of p is greater than the stored cell
            // with the same key. Set j to be the index of u with the lower
            // score.
            uint32_t currScore = iter->second;
            if((int)currScore >= p->G)
            {
                j = i;
            }
            else
            {
                j = iter->second >> 32;
                iter->second = i << 32 | p->G;
            }
        }
        else
        {
            // initialize hash for cell i
            hash[key] = i << 32 | p->G;
        }

        if(j >= 0)
        {
            printf("marking p [%zu %zu] as deleted cidx: %d [DUP]\n", p->interval.lower, p->interval.upper, p->parent_cidx);
            p = &u->cells[j];
            p->interval.lower = -1;
            p->interval.upper = -1;
            p->G = 0;
            if(p->parent_idx >= 0)
                u->cells[p->parent_idx].children_idx[p->parent_cidx] = -3;
            n += 1;
        }
    }
    printf("removed %d duplicate entries\n", n);
}

// Fill the values of C[0] depending on the values in the other 3 cells
int fillCells(const LRParams& params, int match_score, LRCell* c[4])
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
