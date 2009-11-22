//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BWTTraverse - traverse a bwt, outputting all
// strings that are length greater than some threshold
//
#include "BWTTraverse.h"

void BWTTraverse::extract(const BWT* pBWT, unsigned int len)
{
	// Keep a stack of the elements to visit and 
	// a string which is the reverse string corresponding
	// to the current stack elements
	TraverseStack stack;
	std::string rev_str;
	size_t iter_count = 0;
	size_t output_count = 0;

	for(size_t i = 0; i < DNA_ALPHABET_SIZE; ++i)
	{
		char base_char = ALPHABET[i];

		// Initialize the stack and rev_str
		BWTInterval range;
		BWTAlgorithms::initInterval(range, base_char, pBWT);
		AlphaCount ext = BWTAlgorithms::getExtCount(range, pBWT);
		stack.push(TraverseElem(range, ext));
		rev_str.append(1, base_char);

		while(!stack.empty())
		{
			++iter_count;
			TraverseElem& curr = stack.top();
			curr.goNext();
			bool doPop = !curr.isValid(); 

			// If rev_str is long enough output it and stop the traversal at this depth
			if(rev_str.length() == len)
			{
				if(output_count % 1000000 == 0)
					std::cerr << "output: " << output_count << "\n";

				doPop = true;
				size_t multiplicity = curr.getRange().size();
				printf(">%zu %u %zu\n%s\n", output_count, len, multiplicity, reverse(rev_str).c_str());
				++output_count;
			}

			if(doPop)
			{
				// Remove one char from the end of the string 
				// and one element from the stack
				rev_str.erase(rev_str.length() - 1, 1);
				stack.pop();
			}
			else
			{
				// Traverse one level deeper 
				char b = curr.getCurrChar();
				BWTInterval subrange = curr.getRange();
				BWTAlgorithms::updateInterval(subrange, b, pBWT);
				AlphaCount subext = BWTAlgorithms::getExtCount(subrange, pBWT);
				stack.push(TraverseElem(subrange, subext));
				rev_str.append(1, b);
			}
		}
	}

	std::cerr << "Num output: " << output_count << " num iterations: " << iter_count << " count/len: " << (double)iter_count / len << "\n";
}
