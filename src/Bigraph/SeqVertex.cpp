#include "SeqVertex.h"
#include "Util.h"

void SeqVertex::merge(const Vertex* pV2, const Edge& e)
{
	// run the parent merge
	Vertex::merge(pV2, e);

	SeqVertex* pSV2 = (SeqVertex*)pV2;

	// merge the data
	Sequence seq2 = (e.getComp() == EC_SAME) ? pSV2->getSeq() : reverseComplement(pSV2->getSeq());
	Sequence leftSeq;
	Sequence rightSeq;

	if(e.getDir() == ED_SENSE)
	{
		leftSeq = getSeq();
		rightSeq = seq2;
	}
	else
	{
		leftSeq = seq2;
		rightSeq = getSeq();
	}

	// Ensure the overlap is correct
	Sequence leftOverlap = leftSeq.substr(leftSeq.length() - e.getOverlap());
	Sequence rightOverlap = rightSeq.substr(0, e.getOverlap());
	if(leftOverlap != rightOverlap)
	{
		std::cerr << "SeqVertex::merge: left overlap != right overlap (" << leftOverlap << "," << rightOverlap << ")" << std::endl;
		std::cerr << "Edge: " << e << std::endl;
		assert(false && "Merge failure");
	}

	// Get the substring of the right sequence to add to the left
	Sequence rightOverhang = rightSeq.substr(e.getOverlap());
	Sequence merged = leftSeq + rightOverhang;
	assert(merged.length() == leftSeq.length() + rightSeq.length() - e.getOverlap());
	std::cout << "Merged: " << merged << std::endl;
	m_sequence = merged;
}

