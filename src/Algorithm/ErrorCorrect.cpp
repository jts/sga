//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ErrorCorrect - Methods to identify and correct read errors
//
#include "ErrorCorrect.h"

// Perform error correction on the given vertex
std::string ErrorCorrect::correctVertex(Vertex* pVertex, size_t simpleCutoff)
{
	if(pVertex->countEdges() < simpleCutoff)
		return simpleCorrect(pVertex);
	else
		return trieCorrect(pVertex);
}

// simpleCorrect calls a straightforward consensus sequence
// based on all the reads overlapping pVertex
std::string ErrorCorrect::simpleCorrect(Vertex* pVertex)
{
	MultiOverlap mo = pVertex->getMultiOverlap();
	return mo.calculateConsensusFromPartition(0.01);
}

// trieCorrect builds tries from the overlapping reads
// to attempt to account for overcollapsed repeats 
std::string ErrorCorrect::trieCorrect(Vertex* pVertex)
{
	std::string consensus;
	SeqTrie leftTrie;
	SeqTrie rightTrie;
	pVertex->fillTries(0.01, &leftTrie, &rightTrie);
	
	// Re-map low quality branches in the tries
	leftTrie.remodel(2, log(0.01));
	rightTrie.remodel(2, log(0.01));

	std::string original = pVertex->getSeq();

	PathScoreVector left_psv;
	PathScoreVector right_psv;

	leftTrie.score(reverse(original), 0.01, left_psv);
	rightTrie.score(original, 0.01, right_psv);

	size_t leftBestIdx = 0;
	size_t rightBestIdx = 0;
	double leftBestScore = -std::numeric_limits<double>::max();
	double rightBestScore = -std::numeric_limits<double>::max();

	std::cout << "\nLEFTPSV:\n";
	for(size_t i = 0; i < left_psv.size(); ++i)
	{
		left_psv[i].reverse();
		left_psv[i].print();

		std::string path_str = left_psv[i].path_corrected;
		std::string ods = getDiffString(path_str, original);
		printf("CDO: %s\n", ods.c_str());

		if(left_psv[i].path_score > leftBestScore)
		{
			leftBestScore = left_psv[i].path_score;
			leftBestIdx = i;
		}
	}

	std::cout << "\nRIGHTPSV:\n";
	for(size_t i = 0; i < right_psv.size(); ++i)
	{
		right_psv[i].print();
		std::string path_str = right_psv[i].path_corrected;

		std::string ods = getDiffString(path_str, original);
		printf("CDO: %s\n", ods.c_str());

		if(right_psv[i].path_score > rightBestScore)
		{
			rightBestScore = right_psv[i].path_score;
			rightBestIdx = i;
		}
	}
	double cutoff = -40;
	if(left_psv.empty() && right_psv.empty() || (leftBestScore < cutoff && rightBestScore < cutoff))
	{
		consensus = original;
	}
	else if(left_psv.empty() || leftBestScore < cutoff)
	{
		consensus = right_psv[rightBestIdx].path_corrected;
	}
	else if(right_psv.empty() || rightBestScore < cutoff)
	{
		consensus = left_psv[leftBestIdx].path_corrected;
	}
	else
	{
		consensus.reserve(original.size());
		std::cout << "\nComputing consensus from paths\n";
		PathScore& left = left_psv[leftBestIdx];
		PathScore& right = right_psv[rightBestIdx];
		left.print();
		right.print();

		// Compute the combined consensus
		for(size_t i = 0; i < original.size(); ++i)
		{
			if(left.path_corrected[i] == right.path_corrected[i])
			{
				consensus.push_back(left.path_corrected[i]);
			}
			else
			{
				double ls = log(1 - left.path_score) + left.probVector[i];
				double rs = log(1 - right.path_score)+ right.probVector[i];

				if(ls < rs)
					consensus.push_back(left.path_corrected[i]);
				else
					consensus.push_back(right.path_corrected[i]);
			}
		}
	}
	return consensus;
}
