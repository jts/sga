//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// StringEdge class, derived from Bigraph Edge
//
#include "StringEdge.h"
#include "StringVertex.h"

boost::object_pool<StringEdge>* StringEdge::m_spMempool; 

// Compare string edge points by length
bool StringEdgeLenComp::operator()(const Edge* pA, const Edge* pB)
{
	const StringEdge* pSA = static_cast<const StringEdge*>(pA);
	const StringEdge* pSB = static_cast<const StringEdge*>(pB);
	return pSA->getSeqLen() < pSB->getSeqLen();
}

// Flip the string of an edge
void StringEdge::flip()
{
	Edge::flip();
}

// Get the edge's label
std::string StringEdge::getLabel() const
{
	const StringEdge* pTwin = CSE_CAST(getTwin());
	const StringVertex* pEndpoint = CSV_CAST(m_pEnd);
	
	// get the unmatched coordinates in V2
	SeqCoord unmatched = pTwin->getMatchCoord().complement();
	std::string seq = unmatched.getSubstring(pEndpoint->getSeq());

	if(getComp() == EC_REVERSE)
		seq = reverseComplement(seq);

	return seq;
}

// Get the matching portion of V1 described by this edge
std::string StringEdge::getMatchStr() const
{
	StringVertex* pVert = SV_CAST(getStart());
	return m_matchCoord.getSubstring(pVert->getSeq());
}

// Return the length of the sequence
size_t StringEdge::getSeqLen() const
{
	const StringEdge* pTwin = CSE_CAST(getTwin());
	SeqCoord unmatched = pTwin->getMatchCoord().complement();
	return unmatched.length();
}

// return the mapping from V1 to V2 via this edge
// If necessary, flip the SC of V2 so that it is
// in the same coordinate system as V1
Match StringEdge::getMatch() const
{
	const StringEdge* pTwin = CSE_CAST(getTwin());

	const SeqCoord& sc = getMatchCoord();
	const SeqCoord& tsc = pTwin->getMatchCoord();

	return Match(sc, tsc, getComp() == EC_REVERSE, m_numDiff);
}

void StringEdge::validate() const
{
	const StringEdge* pTwin = CSE_CAST(getTwin());
	std::string m_v1 = getMatchStr();
	std::string m_v2 = pTwin->getMatchStr();

	if(getComp() == EC_REVERSE)
		m_v2 = reverseComplement(m_v2);

	bool error = false;
	if(m_v1.length() != m_v2.length())
	{
		std::cerr << "Error, matching strings are not the same length\n";
		error = true;
	}
	else
	{
		int numDiff = countDifferences(m_v1, m_v2, m_v1.length());
		if(numDiff != m_numDiff)
		{
			std::cerr << "Error, number of differences between m1 and m2 does not match expected (" 
			          << numDiff << " != " << m_numDiff << ")\n";
			error = true;
		}
	}

	if(error)
	{
		std::cerr << "V1M: " << m_v1 << "\n";
		std::cerr << "V2M: " << m_v2 << "\n";
		std::cerr << "V1MC: " << getMatchCoord() << "\n";
		std::cerr << "V2MC: " << pTwin->getMatchCoord() << "\n";
		std::cerr << "V1: " << SV_CAST(getStart())->getSeq() << "\n";
		std::cerr << "Validation failed for edge " << *this << "\n";
		assert(false);
	}
}

// Join pEdge into the start of this edge
// This merges the intervals
void StringEdge::join(const Edge* pEdge)
{
	const StringEdge* pSE = CSE_CAST(pEdge);
	//std::cout << "J Edge : " << *pEdge << "\n";
	//std::cout << "J Edge2: " << *this << "\n";
	
	Match m12 = pSE->getMatch();
	Match m23 = getMatch();
	m_matchCoord = m12.inverseTranslate(m23.coord[0]);
	
	// Update the base class members
	Edge::join(pEdge);
}

// Join pEdge into the end of this edge
// This involves merging the intervals
void StringEdge::extend(const Edge* pEdge)
{
	// Update the base class members
	Edge::extend(pEdge);
}

// Update is called when the edge has been modified in some way
// This means its been merged/moved and its difference count should be updated
void StringEdge::update()
{
	updateDifferenceCount();
	Edge::update();
}

void StringEdge::offsetMatch(int offset)
{
	m_matchCoord.interval.start += offset;
	m_matchCoord.interval.end += offset;
}

void StringEdge::extendMatch(int ext_len)
{
	m_matchCoord.interval.end += ext_len;
}

// Bump the edges of the match outwards so it covers the entire 
// sequence
void StringEdge::extendMatchFullLength()
{
	if(m_matchCoord.isLeftExtreme())
		m_matchCoord.interval.end = m_matchCoord.seqlen - 1;
	else
		m_matchCoord.interval.start = 0;
}

void StringEdge::updateSeqLen(int newLen)
{
	m_matchCoord.seqlen = newLen;
}

void StringEdge::updateDifferenceCount()
{
	const StringEdge* pTwin = CSE_CAST(getTwin());
	std::string m_v1 = getMatchStr();
	std::string m_v2 = pTwin->getMatchStr();

	if(getComp() == EC_REVERSE)
		m_v2 = reverseComplement(m_v2);

	m_numDiff = countDifferences(m_v1, m_v2, m_v1.length());
}
