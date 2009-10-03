#ifndef GRAPHCOMMON_H
#define GRAPHCOMMON_H

#include <vector>
#include "Util.h"

// Flags specifying how the dot file should be drawn
enum DotFlags
{
	DF_UNDIRECTED = 0x01,
	DF_ANNOTATIONS = 0x02
};


// Typedefs
typedef std::string VertexID;
typedef std::vector<VertexID> VertexIDVec;


// An edge description is the triplet of values
// that is needed to uniquely identify an edge
struct EdgeDesc
{
	EdgeDesc(VertexID i, EdgeDir d, EdgeComp c) : id(i), dir(d), comp(c) {}
	VertexID id;
	EdgeDir dir;
	EdgeComp comp;

	// Operators
	bool operator<(const EdgeDesc& obj) const
	{
		if(id < obj.id)
			return true;
		else if(id > obj.id)
			return false;
		else if(dir < obj.dir)
			return true;
		else if(dir > obj.dir)
			return false;
		else if(comp < obj.comp)
			return true;
		else if(comp > obj.comp)
			return false;
		return false;
	}

	bool operator==(const EdgeDesc& obj) const
	{
		return id == obj.id && dir == obj.dir && comp == obj.comp;
	}

	friend std::ostream& operator<<(std::ostream& out, const EdgeDesc& ed)
	{
		out << ed.id << "," << ed.dir << "," << ed.comp;
		return out;
	}
};

#endif 
