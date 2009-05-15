#include "UniEst.h"
#include <math.h>

int main(int argc, char** argv)
{
	parseOptions(argc, argv);
	
	// Filenames
	std::string contigFile(argv[optind++]);
	std::string alignFile(argv[optind++]);
	std::string pairedFile(argv[optind++]);

	// Coverage map
	CoverageMap covMap;

	// Read the contigs to create the coverage map
	std::ifstream contigReader(contigFile.c_str());
	assert(contigReader.is_open());
	Contig c;
	while(contigReader >> c)
	{
		UniData ud;
		covMap[c.id] = ud;
		CoverageMap::iterator ret = covMap.find(c.id);
		ret->second.kmerCoverage.resize(c.length - opt::k + 1, 0);
	}
	contigReader.close();
/*
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
			addCoverage(covMap, ka);
		}
	}
	alignReader.close();
*/
	// Read the alignments to generate the coverage
	std::ifstream pairedReader(pairedFile.c_str());
	assert(pairedReader.is_open());
	
	std::string line;
	while(getline(pairedReader, line))
	{
		std::stringstream ss(line);
		std::string read1, read2;
		KAlignment ka1, ka2;
		ss >> read1 >> ka1 >> read2 >> ka2;
		addOverhangCoverage(covMap, ka1, ka2);
	}
	pairedReader.close();

	estimateUniquenessByOverhang(covMap);

/*
	double mean_parameter = estimateGlobalMeanCoverage(covMap);
	std::cerr << "Mean coverage estimated to be: " << mean_parameter << std::endl;

	estimateUniquenessByDepth(covMap, mean_parameter);
*/
}

void estimateUniquenessByOverhang(CoverageMap& covMap)
{
	for(CoverageMap::iterator iter = covMap.begin(); iter != covMap.end(); ++iter)
	{
		std::cout << iter->first << "\t" << iter->second.kmerCoverage.size();
		for(size_t idx = 0; idx <= 1; ++idx)
		{
			StringSet overhang = iter->second.pairedCoverage[idx];
			std::cout << "\t" << overhang.size();
		}
		std::cout << "\n";
	}
}

void estimateUniquenessByDepth(CoverageMap& covMap, double mean_param)
{
	int high_cn = 50;

	for(CoverageMap::iterator iter = covMap.begin(); iter != covMap.end(); ++iter)
	{
		DoubleVec lpVec;
		IntVec& iVec = iter->second.kmerCoverage;
		size_t len = iVec.size();

		for(int cn = 1; cn < high_cn; ++cn)
		{
			double log_p = 0;
			for(size_t i = 0; i < iVec.size(); ++i)
			{
				log_p += log_poisson(iVec[i], cn*mean_param);
			}
			lpVec.push_back(log_p);
		}


		// Select the max value from the log-probability vector
		size_t best_cn = lpVec.size();
		double best_p = 0.0f;
		for(size_t i = 0; i < lpVec.size(); ++i)
		{
			if(best_cn == lpVec.size() || lpVec[i] > best_p)
			{
				best_cn = i + 1;
			}
		}

		sort(lpVec.begin(), lpVec.end());
		double best = lpVec[lpVec.size() - 1];
		double second = lpVec[lpVec.size() - 2];

		std::string ustr;
		if(len <= 25)
		{
			ustr = "TOOSMALL";
		}
		else 
		{
			if(best_cn == 1)
				ustr = "UNIQUE";
			else
				ustr = "REPEAT";
		}

		std::cout << iter->first << "\t" << len << "\t" << best << "\t" << second << "\t" 
				<< ustr << "\t" << best_cn << std::endl;
	}
}

double log_poisson(int k, double m)
{

	double f_k = log_factorial(k);
	double p = (double)k * log(m) - m - f_k;
	//std::cout << "k: " << k << " f: " << f_k << " m: " << m << " p: " << p << std::endl;
	return p;
}

double log_factorial(int k)
{
	double result = 0;

	while(k > 0)
		result += log(k--); //slow
	return result;
}

double estimateGlobalMeanCoverage(CoverageMap& covMap)
{
	std::vector<LenMean> lmVec;
	for(CoverageMap::iterator iter = covMap.begin(); iter != covMap.end(); ++iter)
	{
		// Calculate the mean coverage for this contig
		IntVec& covVec = iter->second.kmerCoverage;
		double curr_mean = getMeanCoverage(covVec);
		LenMean v;
		v.length = covVec.size();
		v.mean = curr_mean;
		lmVec.push_back(v);
	}

	std::sort(lmVec.begin(), lmVec.end(), LenMean::sortByLen);

	std::vector<LenMean> trimmed;
	
	// Pick the median of the top 10%
	size_t low_idx = (size_t)(0.9*(double)lmVec.size());

	for(size_t i = low_idx; i < lmVec.size(); ++i)
	{
		trimmed.push_back(lmVec[i]);
	}

	std::sort(trimmed.begin(), trimmed.end(), LenMean::sortByMean);

	double globalMean;
	size_t h_idx = trimmed.size() >> 1; 
	
	if(trimmed.size() % 2 == 0)
	{
		globalMean = (trimmed[h_idx].mean + trimmed[h_idx-1].mean) / 2;
	}
	else
	{
		globalMean = trimmed[h_idx].mean;
	}
	return globalMean;
}

double getMeanCoverage(const IntVec& iVec)
{
	double sum = 0.f;
	for(size_t i = 0; i < iVec.size(); ++i)
		sum += iVec[i];
	return sum / iVec.size();
}

void printCoverage(CoverageMap& covMap, ContigID id)
{
	CoverageMap::iterator iter = covMap.find(id);
	assert(iter != covMap.end());
	IntVec& covVec = iter->second.kmerCoverage;
	for(size_t i = 0; i < covVec.size(); ++i)
		std::cout << i << " " << covVec[i] << "\n";
}

void addCoverage(CoverageMap& covMap, const KAlignment ka)
{
	CoverageMap::iterator iter = covMap.find(ka.contig_id);
	assert(iter != covMap.end());
	IntVec& covVec = iter->second.kmerCoverage;

	size_t start = ka.contig_start_pos;
	size_t end = start + (ka.align_length - opt::k + 1);
	assert(start < covVec.size() && end <= covVec.size());
	for(size_t i = start; i < end; ++i)
		covVec[i]++;
}

void addOverhangCoverage(CoverageMap& covMap, const KAlignment& ka1, const KAlignment& ka2)
{
	CoverageMap::iterator iter = covMap.find(ka1.contig_id);
	assert(iter != covMap.end());

	size_t start = ka2.contig_start_pos;
	size_t end = start + (ka2.align_length - opt::k + 1);
	size_t idx = ka1.is_reverse;

	for(size_t i = start; i < end; ++i)
	{
		std::stringstream ustr;
		ustr << ka2.contig_id << ":" << i;
		iter->second.pairedCoverage[idx].insert(ustr.str());
		//std::cout << "posstr: " << ustr.str() << std::endl;
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

