#ifndef GRAPHCOMMON_H
#define GRAPHCOMMON_H

#include <vector>

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



#endif 
