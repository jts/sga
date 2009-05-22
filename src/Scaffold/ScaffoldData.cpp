#include "ScaffoldData.h"
#include <iterator>

void ScaffoldData::addLink(SLink sl, size_t idx)
{
	m_linkVecs[idx].push_back(sl);
}

std::ostream& operator<<(std::ostream& out, const ScaffoldData& sd)
{
	out << sd.m_contig.getID() << " LINKS: \n";
	for(size_t idx = 0; idx <= 1; ++idx)
	{
		const SLinkVec& lv = sd.m_linkVecs[idx];
		out << "\t";
		std::copy(lv.begin(), lv.end(), std::ostream_iterator<SLink>(std::cout, "\t"));
		out << "\n";
	}
	return out;
}

