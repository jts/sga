//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// mkqs - multikey quicksort
//
// Perform a ternary quicksort of strings as described in
// Bentley and Sedgewick, 1997
//
// Example code was downloaded from http://www.cs.princeton.edu/~rs/strings/demo.c
// Modified by JTS to take in a comparator and use a generic type
//

#ifndef MKQS_H
#define MKQS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#define mkqs_swap(a, b) { T tmp = x[a]; x[a] = x[b]; x[b] = tmp; }

// Swap [i..i+n] and [j..j+n] in x
template<typename T>
inline void vecswap(int i, int j, int n, T* x)
{   
	while (n-- > 0) 
	{
		mkqs_swap(i, j);
        i++;
        j++;
    }
}

template<typename T, typename PrimarySorter, typename FinalSorter>
void mkqs(T* x, int n, int depth, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{   
	int a, b, c, d, r, v;
    if (n <= 1)
        return;
	
	// Select a random element to use as the key
	// and swap it to the front
    a = rand() % n;
    mkqs_swap(0, a);
	// Return the character at position depth
	// which will be used for the split
	v = primarySorter.getChar(x[0],depth);
	
	// Swap all the elements that have a character at pos depth
	// that is less than v to the low portion of the array
	// and all that are higher to the upper portion
	a = b = 1;
    c = d = n-1;
    for (;;) 
	{
        while (b <= c && (r = primarySorter.getChar(x[b],depth) - v) <= 0)
		{
            if (r == 0) 
			{
				mkqs_swap(a, b); 
				a++; 
			}
            b++;
        }
        while (b <= c && (r = primarySorter.getChar(x[c],depth) - v) >= 0) 
		{
            if (r == 0) 
			{ 
				mkqs_swap(c, d); 
				d--; 
			}
            c--;
        }
        if (b > c) 
			break;
        mkqs_swap(b, c);
        b++;
        c--;
    }


    r = std::min(a, b-a);
	vecswap(0, b-r, r, x);

    r = std::min(d-c, n-d-1); 
	vecswap(b, n-r, r, x);

    r = b-a; 

	mkqs(x, r, depth, primarySorter, finalSorter);

    if (primarySorter.getChar(x[r], depth) != 0)
	{
        mkqs(x + r, a + n-d-1, depth+1, primarySorter, finalSorter);
	}
	else
	{
		// Finalize the sort by using std::sort
		int n2 = a + n - d - 1;
		std::sort(x + r, x + r + n2, finalSorter);
	}

    r = d-c; 
	mkqs(x + n-r, r, depth, primarySorter, finalSorter);
}

#define mkqs_swap2(a, b) { T t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(p) (primarySorter.getChar(*(p), depth))
#define elem2char(e, d) (primarySorter.getChar((e), (d)))

template<typename T>
inline void vecswap2(T* a, T* b, int n)
{   while (n-- > 0) 
	{
        T t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

template<typename T, typename PrimarySorter>
T* med3func(T* a, T* b, T* c, int depth, const PrimarySorter& primarySorter)
{   int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3(a, b, c) med3func(a, b, c, depth, primarySorter)
template<typename T, typename PrimarySorter, typename FinalSorter>
inline void inssort(T* a, int n, int d, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{   
	T *pi, *pj, s, t;
    for (pi = a + 1; --n > 0; pi++)
	{
        for (pj = pi; pj > a; pj--) 
		{
			/*
            // Inline strcmp: break if *(pj-1) <= *pj
			T s = *(pj - 1);
			T t = *pj;
			char* 
			int i = 0;
			char curr_s;
			char curr_t;
			while(true)
			{
				curr_s = elem2char(s, d + i);
				curr_t = elem2char(t, d + i);
				if(curr_s != curr_t || curr_s == 0)
					break;
				++i;
			}
            if (curr_s < curr_t || (curr_s == curr_t && finalSorter(s, t))) 
                break;
            mkqs_swap2(pj, pj-1);
			*/
            // Inline strcmp: break if *(pj-1) <= *pj
			T elem_s = *(pj - 1);
			T elem_t = *pj;
			const char* s = primarySorter.getChrPtr(elem_s);
			const char* t = primarySorter.getChrPtr(elem_t);

			for (s=s+d, t=t+d; *s==*t && *s!=0; s++, t++)
                ;
            if (*s < *t || (*s == *t && finalSorter(elem_s, elem_t)))
                break;
            mkqs_swap2(pj, pj-1);
    	}
	}
}


template<typename T, typename PrimarySorter, typename FinalSorter>
void mkqs2(T* a, int n, int depth, const PrimarySorter& primarySorter, const FinalSorter& finalSorter)
{   
	int r, partval;
    T *pa, *pb, *pc, *pd, *pl, *pm, *pn, t;
   
   	
   	if (n < 10) 
	{
        inssort(a, n, depth, primarySorter, finalSorter);
        return;
    }
	

    pl = a;
    pm = a + (n/2);
    pn = a + (n-1);

	/*
    if (n > 30) 
	{ 
		// On big arrays, pseudomedian of 9
        d = (n/8);
        pl = med3(pl, pl+d, pl+2*d);
        pm = med3(pm-d, pm, pm+d);
        pn = med3(pn-2*d, pn-d, pn);
    }
	*/
    //pm = med3(pl, pm, pn);
	int mid_idx = rand() % n;

	pm = &a[mid_idx];
    mkqs_swap2(a, pm);
    partval = ptr2char(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;
    for (;;) 
	{
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) 
		{
            if (r == 0) { mkqs_swap2(pa, pb); pa++; }
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) 
		{
            if (r == 0) { mkqs_swap2(pc, pd); pd--; }
            pc--;
        }
        if (pb > pc) break;
        mkqs_swap2(pb, pc);
        pb++;
        pc--;
    }
    pn = a + n;
    r = std::min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
    if ((r = pb-pa) > 1)
        mkqs2(a, r, depth, primarySorter, finalSorter);
    if (ptr2char(a + r) != 0)
        mkqs2(a + r, pa-a + pn-pd-1, depth+1, primarySorter, finalSorter);
	else
	{
		int n2 = pa - a + pn - pd - 1;
		std::sort(a + r, a + r + n2, finalSorter);
	}
    if ((r = pd-pc) > 1)
        mkqs2(a + n-r, r, depth, primarySorter, finalSorter);
}
#endif

