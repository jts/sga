#include "Chain.h"

//
// Statics
//
static const char* IDFIELD = "ID";
static const char* COMPONENTSFIELD = "CM";
static const char* CONTIGIDFIELD = "CI";
static const char* DISTFIELD = "DS";
static const char* ORIENTATIONFIELD = "OR";

//
//
//
void Chain::addLink(ContigID id, int dist, bool isRC)
{
	ChainLink cl(id, dist, isRC);
	m_chain.push_back(cl);
}

//
//
//
std::ostream& operator<<(std::ostream& out, const Chain& c)
{
	std::copy(c.m_chain.begin(), c.m_chain.end(), std::ostream_iterator<ChainLink>(std::cout, " -> "));
	return out;
}

//
//
//
std::ostream& writeCAF(std::ostream& out, const Chain& c)
{
	std::vector<std::string> fields;
	fields.push_back(makeKeyValue(IDFIELD, c.m_id));
	fields.push_back(makeKeyValue(COMPONENTSFIELD, c.m_chain.size()));
	std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(out, "\t"));
	out << "\n";
	for(CLList::const_iterator iter = c.m_chain.begin(); iter != c.m_chain.end(); ++iter)
	{
		out << "\t";
		writeCAF(out, *iter) << "\n";
	}
	return out;
}

//
//
//
std::ostream& writeCAF(std::ostream& out, const ChainLink& c)
{
	std::vector<std::string> fields;
	fields.push_back(makeKeyValue(CONTIGIDFIELD, c.id));
	fields.push_back(makeKeyValue(DISTFIELD, c.dist));
	fields.push_back(makeKeyValue(ORIENTATIONFIELD, c.isRC));
	std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(out, "\t"));
	return out;
}

