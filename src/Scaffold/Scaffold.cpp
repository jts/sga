#include "Scaffold.h"
#include "ScaffoldData.h"
#include "SeqGraph.h"
#include <math.h>
#include <iterator>

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

	buildChains(sdMap);
}


//
//
//
void buildChains(SDMap& sdMap)
{
	// Check if the links are ambigious
	size_t ambi_cutoff = 50;

	ChainVector chainVec;
	SeqGraph scaffoldGraph;

	for(SDMap::iterator sdIter = sdMap.begin(); sdIter != sdMap.end(); ++sdIter)
	{
		ScaffoldData& sd = sdIter->second;
		Contig& sourceContig = sd.getContig();

		if(!sourceContig.isUnique())
		{
			continue;
		}

		// Build the chain if its not ambigious
		Chain chain;

		ContigPositionVector cpVec;

		// Add the root contig
		ContigPosition cpos(sourceContig.getID(), ED_SENSE, EC_SAME, 0, 0 + sourceContig.getLength());
		cpVec.push_back(cpos);
		
		// add the left links
		SLinkVec leftVec = sd.getLinks(0);
		for(SLinkVec::iterator vIter = leftVec.begin(); vIter != leftVec.end(); ++vIter)
		{
			Contig& ctg = getContig(sdMap, vIter->linkedID);
			int end = -vIter->dist;
			int start = end - ctg.getLength();

			EdgeComp ec;
			EdgeDir ed;
			if(vIter->isRC)
			{
				ec = EC_REVERSE;
				ed = ED_ANTISENSE;
			}
			else
			{
				ec = EC_SAME;
				ed = ED_SENSE;
			}

			ContigPosition cpos(vIter->linkedID, ed, ec, start, end);
			cpVec.push_back(cpos);
		}

		// add the right links
		SLinkVec rightVec = sd.getLinks(1);
		for(SLinkVec::iterator vIter = rightVec.begin(); vIter != rightVec.end(); ++vIter)
		{
			Contig& ctg = getContig(sdMap, vIter->linkedID);
			int start = sourceContig.getLength() + vIter->dist;
			int end = start + ctg.getLength();
		
			EdgeComp ec;
			EdgeDir ed;
			if(vIter->isRC)
			{
				ec = EC_REVERSE;
				ed = ED_SENSE;
			}
			else
			{
				ec = EC_SAME;
				ed = ED_ANTISENSE;
			}

			ContigPosition cpos(vIter->linkedID, ed, ec, start, end);
			cpVec.push_back(cpos);
		}

		sort(cpVec.begin(), cpVec.end(), ContigPosition::sortStart);

		// Check whether the layout is unambigious
		bool isAmbigious = false;
		Range prevCoords = cpVec.front().pos;
		for(ContigPositionVector::iterator iter = cpVec.begin() + 1; iter != cpVec.end(); ++iter)
		{
			Range currCoords = iter->pos;
			Range intersection = intersect(prevCoords, currCoords);
			if(intersection.size() > ambi_cutoff)
			{
				isAmbigious = true;
				break;
			}
			prevCoords = currCoords;
		}

		// The contigs are now laid out and ordered, build the scaffold graph
		if(!isAmbigious)
		{
			ContigPositionVector::iterator prev = cpVec.begin();
			addVertexToScaffoldGraph(scaffoldGraph, prev->id);
			for(ContigPositionVector::iterator iter = cpVec.begin() + 1; iter != cpVec.end(); ++iter)
			{
				// We infer the relative complementarity from the comp's wrt to the source node
				EdgeComp relativeComp = (prev->comp != iter->comp) ? EC_REVERSE : EC_SAME;
				
				int dist = iter->pos.start - prev->pos.end;
				addVertexToScaffoldGraph(scaffoldGraph, iter->id);
				addEdgeToScaffoldGraph(scaffoldGraph, prev->id, iter->id, iter->dir, relativeComp, dist);
				
				prev = iter;
			}
		}
		std::string ambiStr = isAmbigious ? "AMBI" : "CLEAN";

		//std::cout << sd;
		//std::cout << "\t\t";
		//std::copy(cpVec.begin(), cpVec.end(), std::ostream_iterator<ContigPosition>(std::cout, "\t"));
		//std::cout << ambiStr << "\n";
	}

	scaffoldGraph.writeDot("scaffold.dot");
}

//
//
//
void addVertexToScaffoldGraph(SeqGraph& graph, ContigID id)
{
	if(!graph.hasVertex(id))
	{
		Vertex* pVertex = new Vertex(id);
		graph.addVertex(pVertex);
	}
}

//
//
//
void addEdgeToScaffoldGraph(SeqGraph& graph, ContigID id1, ContigID id2, EdgeDir dir, EdgeComp comp, int dist)
{
	Edge e(id1, id2, dir, comp, dist);
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

