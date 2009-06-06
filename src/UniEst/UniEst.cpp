#include "UniEst.h"
#include <math.h>

int main(int argc, char** argv)
{
	parseOptions(argc, argv);
	
	// Filenames
	std::string contigFile(argv[optind++]);
	std::string alignFile(argv[optind++]);

	// Coverage map
	CoverageMap covMap;

	// Read the contigs and create the coverage map
	std::ifstream contigReader(contigFile.c_str());
	assert(contigReader.is_open());
	Contig c;
	while(readFasta(contigReader,c))
	{
		UniData ud(c, opt::k);
		covMap.insert(std::make_pair(c.getID(), ud));
	}
	contigReader.close();

	// Read the alignments to generate the coverage
	std::ifstream alignReader(alignFile.c_str());
	assert(alignReader.is_open());
	
	std::string line;
	while(getline(alignReader, line))
	{
		std::stringstream ss(line);
		std::string read_name;
		ss >> read_name;
		KAlignment ka;
		while(ss >> ka)
		{
			CoverageMap::iterator iter = covMap.find(ka.contig_id);
			assert(iter != covMap.end());
			iter->second.addCoverage(opt::k, ka);
		}
	}
	alignReader.close();
	
	double mean_parameter = estimateGlobalMeanCoverage(covMap);
	std::cerr << "Mean coverage estimated to be: " << mean_parameter << std::endl;
	estimateUniquenessByDepth(covMap, mean_parameter);

	// Output contigs
	for(CoverageMap::iterator iter = covMap.begin(); iter != covMap.end(); ++iter)
	{
		UniData& data = iter->second;
		writeCAF(std::cout, data.getContig()) << std::endl;
	}
}

void estimateUniquenessByDepth(CoverageMap& covMap, double mean_param)
{
	int high_cn = 10;
	double threshold = 50;
	for(CoverageMap::iterator iter = covMap.begin(); iter != covMap.end(); ++iter)
	{
		DoubleVec pVec;
		UniData& data = iter->second;
		for(int cn = 1; cn < high_cn; ++cn)
		{
			double p1 = log_poisson(data.getNumReads(), cn * mean_param * data.getContigLen());
			//std::cout << mean_param << " " << data.getContigLen() << " " << data.getNumReads() << "\n";
			//double p = mean_param * data.getContigLen() - data.getNumReads() * log(2);
			pVec.push_back(p1);
		}


		// Select the max value from the log-probability vector
		size_t best_cn = pVec.size();
		double best_p = 0.0f;
		for(size_t i = 0; i < pVec.size(); ++i)
		{
			if(best_cn == pVec.size() || pVec[i] > best_p)
			{
				best_cn = i + 1;
			}
		}

		sort(pVec.begin(), pVec.end());
		double best = pVec[pVec.size() - 1];
		double second = pVec[pVec.size() - 2];

		double ratio = best - second;
		
		
		if(data.getContigLen() <= 50 || ratio < threshold)
		{
			data.getContig().setUniqueFlag(UF_NOCALL);
		}
		else
		{
			if(best_cn == 1)
				data.getContig().setUniqueFlag(UF_UNIQUE);
			else
				data.getContig().setUniqueFlag(UF_REPEAT);
		}

		
		std::cerr << iter->first << "\t" << data.getContigLen() << "\t" << best << "\t" << second << "\t" 
					<< ratio << "\t" << "COPYNUM: " << "\t" << best_cn << std::endl;
		
	}
}

double estimateGlobalMeanCoverage(CoverageMap& covMap)
{
	size_t cutoff = 1000;
	
	int sumReads = 0;
	int sumLen = 0;
	for(CoverageMap::iterator iter = covMap.begin(); iter != covMap.end(); ++iter)
	{
		UniData& data = iter->second;
		if(data.getContigLen() > cutoff)
		{
			sumReads += data.getNumReads();
			sumLen += data.getContigLen();
		}
	}

	return (double)sumReads/(double)sumLen;
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
			case 'k': arg >> opt::k; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (opt::k <= 0) 
	{
		std::cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
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

