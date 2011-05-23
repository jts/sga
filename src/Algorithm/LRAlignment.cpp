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
#include "MultiAlignment.h"
#include "stdaln.h"

//#define BWA_COMPAT_DEBUG 1
//#define BWA_COMPAT_DEBUG_RESOLVE 1

namespace LRAlignment
{
static const int MINUS_INF = -0x3fffffff;
#define MASK_LEVEL 0.90f

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
void bwaswAlignment(const std::string& query, const BWT* pTargetBWT, const SampledSuffixArray* pTargetSSA)
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
    LRHitVector hitsVector(2*query.size());

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

#ifdef BWA_COMPAT_DEBUG
        printf("Outer stack loop [%zu %zu] old_n: %zu\n", old_n, v->interval.lower, v->interval.upper);
#endif
        // TODO: bandwidth test and max depth
        
        // Calculate occurrence values for children of the current entry
        for(int qci = 0; qci < DNA_ALPHABET::size; ++qci)
        {
            char query_child_base = DNA_ALPHABET::getBase(qci);
            BWTInterval child_interval = v->interval;
            BWTAlgorithms::updateInterval(child_interval, query_child_base, pQueryBWT);

#ifdef BWA_COMPAT_DEBUG
            printf("Query child interval(%d): [%zu %zu]\n", qci, child_interval.lower, child_interval.upper);
#endif

            if(!child_interval.isValid())
                continue;

            // Update the hashed count of the number of times this substring needs to be visited
            uint64_t key = child_interval.lower << 32 | child_interval.upper;
            LRHash::iterator hashIter = dawgHash.find(key);
            assert(hashIter != dawgHash.end());
            assert((uint32_t)hashIter->second > 0);
            --hashIter->second;

#ifdef BWA_COMPAT_DEBUG
			printf("hash k1: %zu k2: %zu v1: %zu v2: %d\n", child_interval.lower,  child_interval.upper, hashIter->second >> 32, (uint32_t)hashIter->second);
#endif
            // Create an array of cells for the current node holding
            // scores between the current node and nodes in the prefix trie
            LRStackEntry* u = new LRStackEntry;
            u->interval = child_interval;

            // Loop over the nodes in v
            for(size_t i = 0; i < v->cells.size(); ++i)
            {
                LRCell* p = &v->cells[i];

#ifdef BWA_COMPAT_DEBUG
                printf("processing p = v[%zu] [%zu %zu] is deleted? %d\n", 
                        i, p->interval.lower == -1 ? 0 : p->interval.lower, 
                           p->interval.upper == -1 ? 0 : p->interval.upper, 
                           p->interval.upper == -1);
#endif

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
#ifdef BWA_COMPAT_DEBUG
                    printf("parent has been visited\n");
#endif
                    int upos = v->cells[p->parent_idx].u_idx;

#ifdef BWA_COMPAT_DEBUG
                    printf("    upos: %d\n", upos);
#endif
                    c[1] = upos >= 0 ? &u->cells[upos] : NULL;
                    c[2] = p;
                    c[3] = &v->cells[p->parent_idx];

                    int match_score = qci == p->parent_cidx ? params.match : -params.mismatch;
                    int score = fillCells(params, match_score, c);

#ifdef BWA_COMPAT_DEBUG
                    printf("    score: %d\n", score);
#endif
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
#ifdef BWA_COMPAT_DEBUG
                    printf("parent not visited\n");
#endif
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
#ifdef BWA_COMPAT_DEBUG
                    printf("adding x to u\n");
#endif
                    // set the remaining fields in the current cell
                    x.clearChildren();
                    x.parent_cidx = p->parent_cidx;
                    x.revTargetString = p->revTargetString;
                    x.interval = p->interval;
                    x.q_len = p->q_len + 1;
                    x.t_len = p->t_len;
                    u->cells.push_back(x);

                    // TODO: some heap stuff?
                }

                // Check if we should descend into another node of the prefix trie of the target
#ifdef BWA_COMPAT_DEBUG
                printf("x->G: %d qr: %d i: %zu old_n: %zu\n", x.G, params.gap_open_extend, i, old_n);
#endif
                if( (x.G > params.gap_open_extend /* && heap test */) || i < old_n)
                {
#ifdef BWA_COMPAT_DEBUG
                    printf("    pass test 1\n");
#endif
                    if(p->hasUninitializedChild())
                    {
#ifdef BWA_COMPAT_DEBUG
                        printf("    pass test 2\n");
#endif
                        for(int tci = 0; tci < DNA_ALPHABET::size; ++tci)
                        {
                            if(p->children_idx[tci] != -1)
                            {
#ifdef BWA_COMPAT_DEBUG
                                printf("    child exists at pos %d\n", p->children_idx[tci]);
#endif
                                continue; // already added
                            }
                            else
                            {
#ifdef BWA_COMPAT_DEBUG
                                printf("    child does not exist %d\n", p->children_idx[tci]);
#endif
                            }
                            char target_child_base = DNA_ALPHABET::getBase(tci);
                            BWTInterval target_child_interval = p->interval;
                            BWTAlgorithms::updateInterval(target_child_interval, target_child_base, pTargetBWT);
#ifdef BWA_COMPAT_DEBUG
                            printf("Target child interval(%d)): [%zu %zu]\n", tci, target_child_interval.lower, target_child_interval.upper);
#endif
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

                            y.revTargetString = p->revTargetString + target_child_base;

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
#ifdef BWA_COMPAT_DEBUG
            printf("cnt: %d pos: %d\n", count, position);
#endif
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

#ifdef BWA_COMPAT_DEBUG
                    printf("merging stack entries\n");
#endif
                    mergeStackEntries(w, u);
                }

                if(count == 0)
                {
                    // this node in the dawg will not be visited again
                    // move the stack entry from the pending list to the stack
                    removeDuplicateCells(w, dupHash);
#ifdef BWA_COMPAT_DEBUG
                    printf("moving w from pending to stack[%zu %zu]\n", w->interval.lower, w->interval.upper);
#endif
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
#ifdef BWA_COMPAT_DEBUG
                    printf("saving u to pending [%zu %zu]\n", u->interval.lower, u->interval.upper);
#endif
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
#ifdef BWA_COMPAT_DEBUG
                printf("pushing u to stack [%zu %zu]\n", u->interval.lower, u->interval.upper);
#endif
            }
        } // for qci
        
        // done with v
        delete v;
    } // for all stack

    assert(num_pending == 0);

    //
    resolveDuplicateHits(pTargetBWT, pTargetSSA, hitsVector, 3);

    generateCIGAR(query, params, hitsVector);

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
                q->targetString.assign(p->revTargetString.rbegin(), p->revTargetString.rend());
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
#ifdef BWA_COMPAT_DEBUG
            printf("marking p [%zu %zu] as deleted cidx: %d [DUP]\n", p->interval.lower, p->interval.upper, p->parent_cidx);
#endif
            p = &u->cells[j];
            p->interval.lower = -1;
            p->interval.upper = -1;
            p->G = 0;
            if(p->parent_idx >= 0)
                u->cells[p->parent_idx].children_idx[p->parent_cidx] = -3;
            n += 1;
        }
    }
#ifdef BWA_COMPAT_DEBUG
    printf("removed %d duplicate entries\n", n);
#endif
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

int resolveDuplicateHits(const BWT* pTargetBWT, const SampledSuffixArray* pTargetSSA, LRHitVector& hits, int IS)
{
    if(hits.empty())
        return 0;
    assert(pTargetBWT != NULL);

    if(pTargetBWT != NULL && pTargetSSA != NULL)
    {
        // Convert each hit to target coordinates
        int old_n = hits.size();
        int n = 0;
        
        LRHitVector newHits;
        for(size_t i = 0; i < hits.size(); ++i)
        {
            LRHit* p = &hits[i];
            if(p->interval.isValid() && p->interval.size() <= IS)
                n += p->interval.size();
            else if(!p->interval.isValid()) // bwa compatibility hack
                n += 1; 
            else if(p->G > 0)
                ++n;
        }
#ifdef BWA_COMPAT_DEBUG_RESOLVE
        printf("Total hits: %zu\n", hits.size());
        printf("Reallocating hits array to size: %d\n", n);
#endif
        newHits.resize(n);

        int j = 0;
        for(int i = 0; i < old_n; ++i)
        {
            LRHit* p = &hits[i];
            if(p->interval.isValid() && p->interval.size() <= IS)
            {
                for(int64_t k = p->interval.lower; k <= p->interval.upper; ++k)
                {
                    newHits[j] = *p;
                    SAElem elem = pTargetSSA->calcSA(k, pTargetBWT);
                    newHits[j].targetID = elem.getID();
                    newHits[j].position = elem.getPos();
                    newHits[j].interval.lower = 0;
                    newHits[j].interval.upper = -1;

#ifdef BWA_COMPAT_DEBUG_RESOLVE
                    printf("Created new hit at position %d\n", newHits[j].position);
#endif
                    ++j;
                }
            }
            else if(p->G > 0)
            {
                newHits[j] = *p;
                SAElem elem = pTargetSSA->calcSA(p->interval.lower, pTargetBWT);
                newHits[j].targetID = elem.getID();
                newHits[j].position = elem.getPos();
                newHits[j].interval.lower = 0;
                newHits[j].interval.upper = -1;
                newHits[j].flag |= 1;
#ifdef BWA_COMPAT_DEBUG_RESOLVE
                printf("Created new hit at position %d\n", newHits[j].position);
#endif
                ++j;
            }
        }

        // Swap the new hits structure with the old hits
        hits.swap(newHits);
    }

    // sort hits by score
    std::sort(hits.begin(), hits.end(), LRHit::compareG);
    
    int i,j;
    // Remove hits with significant overlaps on the query or target sequences
    for(i = 1; i < (int)hits.size(); ++i)
    {
        LRHit* p = &hits[i];
        if(p->G == 0)
            break;
        for(j = 0; j < i; ++j)
        {
            LRHit* q = &hits[j];
            bool compatible = true;

            // hit q has already been removed
            if(q->G == 0)
                continue;

            // if the hits are to different targets they can't overlap
            if(p->targetID != q->targetID)
                continue;

            if(p->interval.upper == -1 && q->interval.upper == -1)
            {
                // Calculate the overlap length on the query
                int qol = (p->end < q->end? p->end : q->end) - (p->beg > q->beg? p->beg : q->beg);
                if(qol < 0)
                    qol = 0;

                if((float)qol / (p->end - p->beg) > MASK_LEVEL || (float)qol / (q->end - q->beg) > MASK_LEVEL) 
                {
                    int64_t tol = 0;
                    if(p->position + p->length < q->position + q->length)
                        tol = p->position + p->length;
                    else
                        tol = q->position + q->length;

                    if(p->position > q->position)
                        tol -= p->position;
                    else
                        tol -= q->position;
                    
					if((double)tol / p->length > MASK_LEVEL || (double)tol / q->length > MASK_LEVEL)
						compatible = false;

#ifdef BWA_COMPAT_DEBUG_RESOLVE
                    printf("    idx = (%d,%d) G=(%d,%d) qol: %d tol: %zu\n", i,j, p->G, q->G, qol, tol);
#endif
                    if(!compatible)
                    {
                        p->G = 0;
                        break;
                    }
                }
            }
        }
    }

    int new_n = i;

    for(i = j = 0; i < new_n; ++i)
    {
        if(hits[i].G == 0)
            continue;
        if(i != j)
            hits[j++] = hits[i];
        else
            ++j;
    }
    hits.resize(j);
#ifdef BWA_COMPAT_DEBUG_RESOLVE
    printf("Final number of hits after duplicate removal: %zu\n", hits.size());
#endif
    return hits.size();
}

void generateCIGAR(const std::string& query, LRParams& params, LRHitVector& hits)
{
    // Set up stdaln data structures
	size_t max_target = ((query.size() + 1) / 2 * params.match + params.gap_ext) / params.gap_ext + query.size(); // maximum possible target length
    path_t* path = (path_t*)calloc(max_target + query.size(), sizeof(path_t));
    AlnParam par;
    int matrix[25];
    par.matrix = matrix;

    // Set matrix
	for (size_t i = 0; i < 25; ++i) par.matrix[i] = -params.mismatch;	
	for (size_t i = 0; i < 4; ++i) par.matrix[i*5+i] = params.match; 
	par.gap_open = params.gap_open;
    par.gap_ext = params.gap_ext;
    par.gap_end = params.gap_ext;
	par.row = 5; 
    par.band_width = params.bandwidth;
    MAlignDataVector mAlignVec;

    for(size_t i = 0; i < hits.size(); ++i)
    {
        std::string querySub = query.substr(hits[i].beg, hits[i].end - hits[i].beg);

        //printf("Alignment DP\n");
        //printf("Target (%d): %s\n", hits[i].position, hits[i].targetString.c_str());
        //printf("Query (%d): %s\n", hits[i].beg, querySub.c_str());
        // Convert strings to 0-3 representation, as needed by the stdaln routine
        // STDALN uses the same representation as DNA_ALPHABET
        int ql = querySub.size();
        uint8_t* pQueryT = new uint8_t[ql];
        
        int tl = hits[i].targetString.size();
        assert(tl <= (int)max_target);

        uint8_t* pTargetT = new uint8_t[hits[i].targetString.size()];
    
        for(int j = 0; j < ql; ++j)
            pQueryT[j] = DNA_ALPHABET::getBaseRank(querySub[j]);

        for(int j = 0; j < tl; ++j)
            pTargetT[j] = DNA_ALPHABET::getBaseRank(hits[i].targetString[j]);

        
        int path_len = 0;
		/*int score =*/ aln_global_core(pTargetT, tl, pQueryT, ql, &par, path, &path_len);
        int cigarLen = 0;
		uint32_t* pCigar = aln_path2cigar32(path, path_len, &cigarLen);
        
        //
        MAlignData maData;
        maData.str = hits[i].targetString;
        maData.position = hits[i].beg;
        
        // Build cigar string
        
        // Add initial padding
        maData.expandedCigar.append(maData.position, 'S');
        
        // Add alignment symbols
        for (int j = 0; j != cigarLen; ++j)
            maData.expandedCigar.append(pCigar[j]>>4, "MID"[pCigar[j]&0xf]);

        mAlignVec.push_back(maData);

        delete [] pQueryT;
        delete [] pTargetT;
        free(pCigar);
    }

    if(!mAlignVec.empty())
        MultiAlignment ma(query, mAlignVec);

    free(path);
}

// Convert a dynamic programming path to a set of padded strings
// Algorithm ported from stdaln.c:aln_stdaln_aux
void path2padded(const std::string& s1, const std::string& s2, std::string& out1, std::string& out2, std::string& outm, path_t* path, int path_len)
{
    path_t* p = path + path_len;

    out1.resize(path_len + 1, 'A');
    out2.resize(path_len + 1, 'A');
    outm.resize(path_len + 1, 'A');

    for (int l = 0; p >= path; --p, ++l) {
        switch (p->ctype) 
        {
            case FROM_M: 
                out1[l] = s1[p->i]; 
                out2[l] = s2[p->j];
                outm[l] = (s1[p->i] == s2[p->j] /*&& s1[p->i] != ap->row*/)? '|' : ' ';
                break;
            case FROM_I: 
                out1[l] = '-'; 
                out2[l] = s2[p->j]; 
                outm[l] = ' '; 
                break;
            case FROM_D: 
                out1[l] = s1[p->i]; 
                out2[l] = '-'; 
                outm[l] = ' '; 
                break;
        }
    }    
}

};
