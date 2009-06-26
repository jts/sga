#ifndef SCAFFOLDDATA_H
#define SCAFFOLDDATA_H

#include "Util.h"
#include "Contig.h"

//
// Structs
//
struct SLink
{
	ContigID linkedID;
	int dist;
	unsigned int numPairs;
	double stdDev;
	bool isRC;
	
	friend std::istream& operator>>(std::istream& in, SLink& sl)
	{
		std::string line;
		in >> line;
		if(line.size () != 0)
		{
			StringVec fields = split(line, ',');
			assert(fields.size() == 5);
			sl.linkedID = fields[0];
			sl.dist = atoi(fields[1].c_str());
			sl.numPairs = atoi(fields[2].c_str());
			std::stringstream fparser(fields[3]);
			fparser >> sl.stdDev;
			sl.isRC = atoi(fields[4].c_str());
		}
		return in;
	}

	friend std::ostream& operator<<(std::ostream& out, const SLink& sl)
	{
		out << sl.linkedID << " " << sl.dist << " " << sl.numPairs << " " << sl.stdDev << " " << sl.isRC;
		return out;
	}		
};


struct ContigPosition
{
	ContigPosition(ContigID cid, EdgeDir d, EdgeComp c, int start, int end) : id(cid), dir(d), comp(c), pos(start, end) { }
	ContigID id;
	EdgeDir dir;
	EdgeComp comp;

	Range pos;

	static bool sortStart(const ContigPosition& c1, const ContigPosition& c2)
	{
		return c1.pos.start < c2.pos.start;
	}

	friend std::ostream& operator<<(std::ostream& out, const ContigPosition& cp)
	{
		out << cp.id << " " << cp.pos;
		return out;
	}

};


//
// Typedefs
//
typedef std::vector<SLink> SLinkVec;

//
//
//
class ScaffoldData
{
	public:
		
		ScaffoldData(Contig& c) : m_contig(c) { }
		void addLink(SLink sl, size_t idx);
		SLinkVec& getLinks(size_t idx) { return m_linkVecs[idx]; }
		Contig& getContig() { return m_contig; }
		friend std::ostream& operator<<(std::ostream& out, const ScaffoldData& sd);

	private:
	
		Contig m_contig;
		SLinkVec m_linkVecs[ED_COUNT];
};

#endif

