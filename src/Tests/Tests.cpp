#include <iostream>
#include <list>

template<typename A>
class TempEdge
{
	A data;
};

template<typename V, typename E>
class TempVert
{
	std::list<TempEdge<E> > m_list;
	V data;
};

template<typename V, typename E>
class TempGraph
{
	typedef TempVert<V,E> VertexType;
	typedef typename std::list<VertexType> VertexList;
	VertexList m_vertices;
};

class OverlapData
{
	int overlap;
};

class SequenceVertex
{	
	
};

int main(int argc, char** argv)
{
	(void)argc;
	(void)argv;

	TempEdge<int> blag;
	TempVert<int, int> vert;
	TempGraph<int, int> graph;

	(void)vert;
	(void)blag;
	return 0;
}


