#ifndef SEQGRAPH_H
#define SEQGRAPH_H

#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include "GraphCommon.h"
#include "Vertex.h"

template<typename VT>
class Bigraph
{

	public:
		//
		// Typedefs
		//
		typedef VT VertexType;
		typedef typename VertexType::EdgeType GraphEdgeType;

		typedef std::map<VertexID, VertexType*> VertexPtrMap;
		typedef typename VertexPtrMap::iterator VertexPtrMapIter;
		typedef typename VertexPtrMap::const_iterator VertexPtrMapConstIter;
		
		typedef bool(*VertexVisitFunction)(Bigraph*, VertexType*);
		typedef std::string(*VertexColorFunction)(typename VertexType::VertexData);

		typedef std::vector<GraphEdgeType> EdgeVec;
		typedef typename EdgeVec::iterator EdgeVecIter;
		typedef EdgeVec Path; // alias
		typedef std::vector<Path> PathVector;
		typedef std::vector<VertexID> VertexIDVec;

	
		Bigraph();
		~Bigraph();

		// Add a vertex
		void addVertex(VertexType* pVert);
		
		// Remove a vertex
		void removeVertex(VertexID id);

		// Check if a vertex exists
		bool hasVertex(VertexID id);

		// Get a vertex
		VertexType* getVertex(VertexID id) const;

		// Add an edge
		void addEdge(const GraphEdgeType& e);

		// Remove an edge
		void removeEdge(const GraphEdgeType& e);

		// Merge vertices
		void mergeVertices(VertexID id1, VertexID id2);

		// Simplify the graph by removing transitive edges
		void simplify();

		// Validate that the graph is sane
		void validate();

		// Flip a given vertex
		void flip(VertexID id);

		// Get the IDs of the vertices that do not branch (both sense/antisense degree <= 1)
		VertexIDVec getNonBranchingVertices() const;

		// Get the linear components of a non-branching graph
		PathVector getLinearComponents();

		// Return all the path of nodes that can be linearally reached from this node
		// The path expands in both directions so the first node in the path is not necessarily the source
		Path constructLinearPath(VertexID id);

		// Print simple summary statistics to stdout
		void stats() const;
		
		// Visit each vertex in the graph and perform the visit function
		bool visit(VertexVisitFunction f);

		// Dump the graph to a dot file
		void writeDot(string filename, int dotFlags = 0, VertexColorFunction colorFunc = &VertexBlackFunction) const;

	private:

		// Set/check the colors for the entire graph
		void setColors(VertexColor c);
		bool checkColors(VertexColor c);

		void followLinear(VertexID id, EdgeDir dir, Path& outPath);

		// Merge two vertices along a specified edge
		void mergeAlongEdge(VertexType* pV1, VertexType* pV2, const GraphEdgeType& edge);

		// Vertex collection
		VertexPtrMap m_vertices;
};

#include "Bigraph.impl"

#endif
