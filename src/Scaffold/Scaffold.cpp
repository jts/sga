#include "Scaffold.h"
#include "ScaffoldData.h"
#include "SeqGraph.h"
#include "SeqVertex.h"
#include <math.h>
#include <iterator>

// Fields
static const char* IDFIELD = "ID";
static const char* COMPONENTSFIELD = "CM";
static const char* CONTIGIDFIELD = "CI";
static const char* DISTFIELD = "DS";
static const char* ORIENTATIONFIELD = "OR";

size_t overlapThreshold = 40;

int main(int argc, char** argv)
{
	parseOptions(argc, argv);
	
	// Filenames
	std::string contigFile(argv[optind++]);
	std::string distanceFile(argv[optind++]);

	SDMap sdMap;

	// Read the contigs
	std::ifstream contigReader(contigFile.c_str());
	assert(contigReader.is_open());

	Contig c;
	while(readCAF(contigReader,c))
	{
		ScaffoldData sd(c);
		sdMap.insert(std::make_pair(c.getID(), sd));
	}
	contigReader.close();

	parseLinks(distanceFile, sdMap);

	buildGraph(sdMap);
}

void buildGraph(SDMap& sdMap)
{
	SeqGraph scaffoldGraph;
	for(SDMap::iterator sdIter = sdMap.begin(); sdIter != sdMap.end(); ++sdIter)
	{
		ScaffoldData& sd = sdIter->second;
		Contig& sourceContig = sd.getContig();
		
		if(!sourceContig.isUnique())
		{
			continue;
		}
	
		addVertexToScaffoldGraph(scaffoldGraph, sdMap, sourceContig.getID());
		
		for(size_t idx = 0; idx <= 1; ++idx)
		{
			EdgeDir currSourceDir = (idx == 0) ? ED_ANTISENSE : ED_SENSE;

			SLinkVec& links = sd.getLinks(idx);
			for(SLinkVec::iterator vIter = links.begin(); vIter != links.end(); ++vIter)
			{
				addVertexToScaffoldGraph(scaffoldGraph, sdMap, vIter->linkedID);
				EdgeComp currComp = (vIter->isRC) ? EC_REVERSE : EC_SAME;
				addEdgeToScaffoldGraph(scaffoldGraph, sourceContig.getID(), vIter->linkedID, currSourceDir, currComp, vIter->dist);
			}
		}
	}

	scaffoldGraph.writeDot("before.dot");
	scaffoldGraph.visit(&cutInconsistent);
	scaffoldGraph.writeDot("cutIncon.dot");
	while(scaffoldGraph.visit(&makeTransitive));
	scaffoldGraph.visit(&cutAmbigious);
	scaffoldGraph.writeDot("after.dot");

	//scaffoldGraph.constructLinearPath("1005");
	PathVector paths = scaffoldGraph.getLinearComponents();
	int id = 0;
	std::ofstream out("scaffold.caf");
	for(PathVector::iterator iter = paths.begin(); iter != paths.end(); ++iter)
	{
		writeScaffold(out, id++, *iter);
	}
	out.close();
}

void writeScaffold(ostream& out, int idNum, const Path& path)
{
	if(path.size() == 0)
	{
		//std::cerr << "Warning zero-element scaffold\n";
		return;
	}

	StringVec headerFields;

	std::stringstream id;
	id << "Scaffold_" << idNum;
	
	headerFields.push_back(makeKeyValue(IDFIELD, id.str()));
	headerFields.push_back(makeKeyValue(COMPONENTSFIELD, path.size() + 1));

	// write the header info
	std::copy(headerFields.begin(), headerFields.end(), std::ostream_iterator<std::string>(out, "\t"));
	out << "\n";
	
	// write the first contig
	Edge firstEdge = path.front();
	writeScaffoldNode(out, firstEdge.getStart(), 0, false);
	out << "\n";
	
	for(Path::const_iterator iter = path.begin(); iter != path.end(); ++iter)
	{
		writeScaffoldNode(out, iter->getEnd(), iter->getOverlap(), iter->getComp() == EC_REVERSE);
		out << "\n";
	}
}

void writeScaffoldNode(ostream& out, VertexID id, int dist, bool orientation)
{
	StringVec fields;
	fields.push_back(makeKeyValue(CONTIGIDFIELD, id));
	fields.push_back(makeKeyValue(DISTFIELD, dist));
	fields.push_back(makeKeyValue(ORIENTATIONFIELD, orientation));
	std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(out, "\t"));
}

bool makeTransitive(SeqGraph* pGraph, Vertex* pVertex)
{
	bool modified = false;
	for(int d = 0; d < ED_COUNT; ++d)
	{
		EdgeDir dir = EDGE_DIRECTIONS[d];
		
		// Get the edges for this direction
		EdgeVec edges = pVertex->getEdges(dir);

		// If there is more than one edge in this direction, check if the edges are consistent
		if(edges.size() > 1)
		{
			LSLVec lslVec;

			for(EdgeVecIter iter = edges.begin(); iter != edges.end(); ++iter)
			{		
				Range r = convertEdgeToRange(pGraph, *iter);
				lslVec.push_back(LinearScaffoldLink(*iter, r));
			}

			sort(lslVec.begin(), lslVec.end(), LinearScaffoldLink::sortStarts);

			//
			// Infer distances between the contigs that are linked to the source contig
			// Check if the distances are all consistent 
			//
			bool isConflicted = false;
			size_t numLinks = lslVec.size();
			assert(numLinks > 0);
			LinearScaffoldLink prev = lslVec[0];
			std::vector<int> distances;

			for(size_t idx = 1; idx < numLinks; ++idx)
			{
				LinearScaffoldLink curr = lslVec[idx];
				int d = curr.range.start - prev.range.end;
				size_t overlap = intersect(prev.range, curr.range).size();
				bool contained = overlap == curr.range.size() || overlap == prev.range.size();
				bool isAmbigious = overlap > overlapThreshold;
				isConflicted = isConflicted || contained || isAmbigious;
				distances.push_back(d);
				prev = curr;
			}

			// Modify the graph by adding the transitive edges if its not conflicting
			// Otherwise remove all the links from this vert in this direction
			if(!isConflicted)
			{
				std::cout << "linearizing off of " << pVertex->getID() << std::endl;
				// Update the graph
				// Remove all links to the source vert except the closest
				// Add the inferred links
				// New links are added between the remaining vertices to make a linear chain
				for(size_t idx = 0; idx < numLinks - 1; ++idx)
				{
					LinearScaffoldLink curr = lslVec[idx];
					LinearScaffoldLink next = lslVec[idx+1];

					int dist = distances[idx];

					// The direction of the inferred edge is the reverse
					// of the edge from the curr node to the source
					EdgeDir cDir = !curr.edge.getTwin().getDir();
					EdgeDir nDir = next.edge.getTwin().getDir();
					
					// The edge comp flag is SAME if the edges have the same
					// orientation wrt the source, OPP otherwise
					EdgeComp cComp = (curr.edge.getComp() == next.edge.getComp()) ? EC_SAME : EC_REVERSE;

					// Create the edge
					Edge inferred(curr.edge.getEnd(), next.edge.getEnd(), cDir, cComp, dist);
					assert(inferred.getTwin().getDir() == nDir); // sanity check

					// Delete the edge from the source to next
					pGraph->removeEdge(next.edge);
					pGraph->removeEdge(next.edge.getTwin());

					// Add the new edge
					Vertex* ns = pGraph->getVertex(inferred.getStart());
					if(!ns->hasEdge(inferred))
					{
						pGraph->addEdge(inferred);
						pGraph->addEdge(inferred.getTwin());
					}
					else
					{
						Edge old = ns->getEdge(inferred);
						std::cout << "Inferred edge already exists " << old << " curr= " << old.getOverlap() << " new= " << inferred.getOverlap() << "\n";
					}

					std::cout << "\tcreated edge " << inferred << std::endl;
				}
			}
			else
			{
				for(EdgeVecIter iter = edges.begin(); iter != edges.end(); ++iter)
				{
					pGraph->removeEdge(*iter);
					pGraph->removeEdge(iter->getTwin());
				}
			}
			modified = true; // the graph has been changed
		}
	}	
	return modified;
}

//
// Check whether the set of edges
//
bool areEdgesConsistent(SeqGraph* pGraph, Vertex* pVertex, EdgeDir dir)
{
	// Get the edges for this direction
	EdgeVec edges = pVertex->getEdges(dir);

	if(edges.size() <= 1)
		return true; //trivial case

	// Convert edges to the ranges that they cover
	LSLVec lslVec;
	for(EdgeVecIter iter = edges.begin(); iter != edges.end(); ++iter)
	{		
		Range r = convertEdgeToRange(pGraph, *iter);
		lslVec.push_back(LinearScaffoldLink(*iter, r));
	}
	
	// Sort by start coordinate
	sort(lslVec.begin(), lslVec.end(), LinearScaffoldLink::sortStarts);

	//
	// Infer distances between the contigs that are linked to the source contig
	// Check if the distances are all consistent 
	//
	bool isConflicted = false;
	size_t numLinks = lslVec.size();
	assert(numLinks > 0);
	LinearScaffoldLink prev = lslVec[0];
	
	std::cout << "Vert: " << pVertex->getID() << " " << dir << std::endl;

	for(size_t idx = 1; idx < numLinks; ++idx)
	{
		LinearScaffoldLink curr = lslVec[idx];
		size_t overlap = intersect(prev.range, curr.range).size();
		bool contained = false;//overlap == curr.range.size() || overlap == prev.range.size();
		bool isAmbigious = overlap > overlapThreshold;
		isConflicted = isConflicted || contained || isAmbigious;
		std::cout << "\t" << prev.range << " " << curr.range << " " << contained << " " << isAmbigious << std::endl;
		prev = curr;
	}
	return !isConflicted;
}

bool cutInconsistent(SeqGraph* pGraph, Vertex* pVertex)
{
	bool modified = false;
	bool isInconsistent = !areEdgesConsistent(pGraph, pVertex, ED_SENSE) || 
							!areEdgesConsistent(pGraph, pVertex, ED_ANTISENSE);

	if(isInconsistent)
	{
		std::cout << "Cutting inconsistent vertex " << pVertex->getID() << std::endl;
		for(int d = 0; d < ED_COUNT; ++d)
		{
			EdgeVec edges = pVertex->getEdges(EDGE_DIRECTIONS[d]);
			for(EdgeVecIter iter = edges.begin(); iter != edges.end(); ++iter)
			{		
				pGraph->removeEdge(*iter);
				pGraph->removeEdge(iter->getTwin());
			}
			modified = true;
		}
	}
	return modified;
}


bool cutAmbigious(SeqGraph* pGraph, Vertex* pVertex)
{
	bool modified = false;
	for(int d = 0; d < ED_COUNT; ++d)
	{
		EdgeDir dir = EDGE_DIRECTIONS[d];
		
		// Get the edges for this direction
		EdgeVec edges = pVertex->getEdges(dir);

		// If there is a single edge in this direction, merge the vertices
		// Don't merge singular self edges though
		if(edges.size() > 1)
		{
			for(EdgeVecIter iter = edges.begin(); iter != edges.end(); ++iter)
			{		
				pGraph->removeEdge(*iter);
				pGraph->removeEdge(iter->getTwin());
			}
			modified = true;
		}
	}
	return modified;
}

//
//
//
Range convertEdgeToRange(const SeqGraph* sg, const Edge& e)
{
	SeqVertex* pVert = (SeqVertex*)sg->getVertex(e.getEnd());
	return Range(e.getOverlap(), e.getOverlap() + pVert->getSeq().length());
}

//
//
//
void addVertexToScaffoldGraph(SeqGraph& graph, SDMap& sdMap, ContigID id)
{
	if(!graph.hasVertex(id))
	{
		SDMap::iterator iter = sdMap.find(id);
		assert(iter != sdMap.end());
		Sequence s = iter->second.getContig().getSequence();
		SeqVertex* pVertex = new SeqVertex(id, s);
		graph.addVertex(pVertex);
	}
}

//
//
//
void addEdgeToScaffoldGraph(SeqGraph& graph, ContigID id1, ContigID id2, EdgeDir dir, EdgeComp comp, int dist)
{
	Edge e(id1, id2, dir, comp, dist);
	//std::cerr << "Adding edge " << e << std::endl;
	graph.addEdge(e);
	graph.addEdge(e.getTwin());
}

//
// 
//
Contig& getContig(SDMap& sdMap, ContigID cid)
{
	SDMap::iterator iter = sdMap.find(cid);
	assert(iter != sdMap.end());
	return iter->second.getContig();
}



//
// Parse a distance file 
//
void parseLinks(std::string filename, SDMap& sdMap)
{
	std::ifstream linkReader(filename.c_str());
	assert(linkReader.is_open());
	
	std::string line;
	while(getline(linkReader, line))
	{
		std::stringstream parser(line);
		ContigID sourceID;
		parser >> sourceID;

		// Get an iterator to the data
		SDMap::iterator iter = sdMap.find(sourceID);
		if(iter == sdMap.end())
		{
			std::cout << "could not find " << sourceID << std::endl;
		}

		assert(iter != sdMap.end());
		ScaffoldData& currSD = iter->second;
	
		// Discard the seperator
		char sep;
		parser >> sep;

		// Read the rest of the line which is the links
		std::string remainder;
		getline(parser, remainder);
		
		// Split the line into two parts, one for each direction
		StringVec dirs = split(remainder, '|');
		assert(dirs.size() == 2);
		for(size_t idx = 0; idx <= 1; ++idx)
		{
			std::stringstream d_parser(dirs[idx]);

			SLink sl;
			while(d_parser >> sl)
			{
				// Check if the link is to a unique contig
				Contig& linkedTo = getContig(sdMap, sl.linkedID);
				(void)linkedTo;
				if(linkedTo.isUnique())
				{
					currSD.addLink(sl, EDGE_DIRECTIONS[idx]);
				}
			}
		}
	}
}

// 
// Handle command line arguments
//
void parseOptions(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 2) 
	{
		std::cerr << PROGRAM ": missing arguments\n";
		die = true;
	} 
	else if (argc - optind > 3) 
	{
		std::cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) 
	{
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
}

