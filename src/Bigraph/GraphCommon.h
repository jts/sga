#ifndef GRAPHCOMMON_H
#define GRAPHCOMMON_H

#include <vector>

// Flags specifying how the dot file should be drawn
enum DotFlags
{
	DF_UNDIRECTED = 0x01,
	DF_ANNOTATIONS = 0x02
};


// Typedefs
typedef std::string VertexID;
typedef std::vector<VertexID> VertexIDVec;

//
// Functions
//
template<typename ET>
std::vector<ET> reversePath(const std::vector<ET>& path)
{
	std::vector<ET> out;
    for(typename std::vector<ET>::const_reverse_iterator iter = path.rbegin(); iter != path.rend(); ++iter)
		out.push_back(iter->getTwin());
	return out;
}

// Default vertex color function, returns black for everything
template<typename D>
std::string VertexBlackFunction(D /*d*/)
{
	return "black";
}



#endif 
