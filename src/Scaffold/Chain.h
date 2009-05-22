#ifndef CHAIN_H
#define CHAIN_H

#include "Util.h"
#include "Contig.h"
#include <list>

//
// Structs
//
struct ChainLink
{
	ChainLink(ContigID cid, int d, bool b) : id(cid), dist(d), isRC(b) {}

	friend std::ostream& operator<<(std::ostream& out, const ChainLink& c)
	{
		out << c.id << "," << c.dist << "," << c.isRC; 
		return out;
	}

	friend std::ostream& writeCAF(std::ostream& out, const ChainLink& c);

	ContigID id;
	int dist;
	bool isRC;

};

//
// Typedefs
//
typedef std::list<ChainLink> CLList;

class Chain
{
	public:

		Chain() {}
		void addLink(ContigID id, int dist, bool isRC);
		void setID(ContigID id) { m_id = id; }

		// Output
		friend std::ostream& writeCAF(std::ostream& out, const Chain& c);
		friend std::ostream& operator<<(std::ostream& out, const Chain& c);


	private:
		ContigID m_id;
		CLList m_chain;
};

typedef std::vector<Chain> ChainVector;

#endif

