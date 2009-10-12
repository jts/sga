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
void vecswap(int i, int j, int n, T* x)
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

#endif

