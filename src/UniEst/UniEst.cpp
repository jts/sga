#include "UniEst.h"
#include "FragmentDistribution.h"
#include "StatsCommon.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>
#include <math.h>

int main(int argc, char** argv)
{
	parseOptions(argc, argv);
	
	// Filenames
	std::string contigFile(argv[optind++]);

	if(opt::verbose >= 1)
	{
		std::cerr << "Contigs file: " << contigFile << "\n";
		std::cerr << "Alignments file: " << (opt::alignFile.empty() ? "none" : opt::alignFile) << "\n";
		std::cerr << "Paired aligns file: " << (opt::pairedFile.empty() ? "none" : opt::pairedFile) << "\n";
		std::cerr << "Threshold: " << opt::threshold << "\n";
		std::cerr << "Len cutoff: " << opt::length_cutoff << "\n";
		std::cerr << "Kmer: " << opt::k << "\n";
		std::cerr << "Writing results to: " << opt::outfile << "\n";
		std::cerr << "Using depth data: " << opt::bUseDepth << "\n";
		std::cerr << "Using paired data: " << opt::bUsePairs << "\n";
	}

	// Parse the contigs file to load the initial data 
	IDUniDataMap udMap;
	parseContigs(contigFile, udMap);

	// TODO: Derive this parameter from data
	double overhang_cutoff = 1.5;
	IntDist fragDist;
	
	// If an alignment file has been specified, read in the file and add depth info
	parseAligns(opt::alignFile, udMap);

	// If a pairs file has been specified, read in the histogram and pair aligns
	if(opt::bUsePairs)
	{
		FragmentDistribution fragSizeCounts;
		fragSizeCounts.readFromFile(opt::histFile);
		fragDist = fragSizeCounts.convertToIntDist(0.99); // Convert to a pdf and trim
		parsePairedAligns(opt::pairedFile, udMap);
	}

	// Initially, fit based on the the top 10% of contigs based on size
	UniDataPVec fitVec = filterByPercent(udMap, 0.9);
	double mean_coverage_parameter, error_rate_parameter;
	size_t numIters = 2; 

	for(size_t i = 0; i < numIters; ++i)
	{
		fitParameters(fitVec, fragDist, opt::bUsePairs, mean_coverage_parameter, error_rate_parameter);
		if(opt::verbose >= 1)
		{
			std::cerr << "Mean estimate: " << mean_coverage_parameter 
					  << " error rate estimate: " << error_rate_parameter << std::endl;
		}
		
		for(IDUniDataMap::iterator iter = udMap.begin(); iter != udMap.end(); ++iter)
		{
			// Don't estimate contigs that are too small
			if(iter->second.getContigLen() < opt::length_cutoff)
			{
				iter->second.getContig().setUniqueFlag(UF_NOCALL);
			}
			else
			{
				bool bAllUnique = true;

				if(opt::bUseDepth)
				{
					UniqueFlag uf = iter->second.estimateByDepth(mean_coverage_parameter, opt::threshold);
					if(uf != UF_UNIQUE)
						bAllUnique = false;
				}

				if(opt::bUsePairs)
				{
					UniqueFlag uf = iter->second.estimateByOverhang(fragDist, 
																			mean_coverage_parameter, 
																			error_rate_parameter, 
																			overhang_cutoff);
					if(uf != UF_UNIQUE)
						bAllUnique = false;
				}

				if(bAllUnique)
				{
					iter->second.getContig().setUniqueFlag(UF_UNIQUE);
				}
				else
				{
					iter->second.getContig().setUniqueFlag(UF_REPEAT);
				}
			}
		}

		// Create the new vector of data points to fit based on the unique contigs
		fitVec = filterByUnique(udMap);
	}
	
	// Print the data
	if(opt::verbose > 1)
	{
		std::cout << "name\tlen\tnum\tdensity\tleft_over_actual\tright_over_actual\t" <<
					 "left_over_expect\tright_over_expect\tdepth_uf\toverhang_uf\n";

		for(IDUniDataMap::iterator iter = udMap.begin(); iter != udMap.end(); ++iter)
		{
			UniData& ud = iter->second;
			ud.printStats(fragDist, mean_coverage_parameter, error_rate_parameter);
		}
	}

	// Output contigs
	std::ofstream outContigs(opt::outfile.c_str());
	for(IDUniDataMap::iterator iter = udMap.begin(); iter != udMap.end(); ++iter)
	{
		UniData& data = iter->second;
		writeCAF(outContigs, data.getContig()) << std::endl;
	}
	outContigs.close();
}

UniDataPVec filterByPercent(const IDUniDataMap& udMap, double percent)
{
	IntVec sizeVec;
	for(IDUniDataMap::const_iterator iter = udMap.begin(); iter != udMap.end(); ++iter)
	{
		sizeVec.push_back(iter->second.getContigLen());
	}

	sort(sizeVec.begin(), sizeVec.end());
	size_t idx = (size_t)((double)sizeVec.size() * percent);
	size_t cutoff = sizeVec[idx];
	
	std::cerr << "Cutoff size is " << cutoff << "\n";
	
	// Filter the data
	UniDataPVec fitVec;
	for(IDUniDataMap::const_iterator iter = udMap.begin(); iter != udMap.end(); ++iter)
	{
		const UniData* pData = &iter->second;
		if(pData->getContigLen() > cutoff)
		{
			fitVec.push_back(pData);
		}
	}
	return fitVec;
}

UniDataPVec filterByUnique(const IDUniDataMap& udMap)
{
	// Filter the data
	UniDataPVec fitVec;
	for(IDUniDataMap::const_iterator iter = udMap.begin(); iter != udMap.end(); ++iter)
	{
		const UniData* pData = &iter->second;
		if(pData->getContig().isUnique())
		{
			fitVec.push_back(pData);
		}
	}
	return fitVec;
}

void fitParameters(const UniDataPVec& fitVec, const IntDist& fragDist, bool fitError, double& mean, double& error_rate)
{
	mean = 0;
	error_rate = 0;
	
	std::cerr << "Fitting using " << fitVec.size() << " data points\n";
	// Fit the mean
	int sumReads = 0;
	int sumLen = 0;
	for(UniDataPVec::const_iterator iter = fitVec.begin(); iter != fitVec.end(); ++iter)
	{
		sumReads += (*iter)->getNumReads();
		sumLen += (*iter)->getContigLen();
	}
	mean = (double)sumReads/(double)sumLen;

	// Compute the error rate parameter with an least-squares fit between the observed - expected overhang coverage
	if(fitError)
	{
		double* pXValues = new double[fitVec.size() * 2];
		double* pYValues = new double[fitVec.size() * 2];

		for(size_t i = 0; i < fitVec.size(); ++i)
		{
			const UniData* pData = fitVec[i];
			size_t len = pData->getContigLen();

			double expected = pData->getExpectedOverhangCoverage(fragDist, mean);

			size_t obs0 = pData->getOverhangCoverage(0);
			size_t obs1 = pData->getOverhangCoverage(1);

			double diff0 = (double)(obs0 - expected);
			double diff1 = (double)(obs1 - expected);

			pXValues[2*i] = len;
			pXValues[2*i+1] = len;

			pYValues[2*i] = diff0;
			pYValues[2*i+1] = diff1;
		}

		// perform a least-squares fit using gsl
		double c0, c1, cov00, cov01, cov11, sumsq;
		gsl_fit_linear(pXValues, 1, pYValues, 1, 2*fitVec.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		delete [] pXValues;
		delete [] pYValues;

		error_rate = c1;
	}
}

//
// Parse the contigs file
//
void parseContigs(std::string file, IDUniDataMap& udMap)
{
	std::ifstream reader(file.c_str());
	if(!reader.is_open())
	{
		std::cerr << "Cannot open " << file << std::endl;
		exit(1);
	}

	Contig c;
	while(readFasta(reader,c))
	{
		UniData ud(c, opt::k);
		udMap.insert(std::make_pair(c.getID(), ud));
	}
	reader.close();
}

//
// Parse an alignment file 
//
void parseAligns(std::string file, IDUniDataMap& udMap)
{
	std::ifstream reader(file.c_str());
	if(!reader.is_open())
	{
		std::cerr << "Cannot open " << file << std::endl;
		exit(1);
	}

	// parse the alignments and generate the kmer coverage
	std::string line;
	while(getline(reader, line))
	{
		std::stringstream ss(line);
		std::string read_name;
		ss >> read_name;
		KAlignment ka;
		while(ss >> ka)
		{
			IDUniDataMap::iterator iter = udMap.find(ka.contig_id);
			assert(iter != udMap.end());
			iter->second.addCoverage(opt::k, ka);
		}
	}
	reader.close();	
}

//
// Parse a paired alignment file
//
void parsePairedAligns(std::string file, IDUniDataMap& udMap)
{
	std::ifstream reader(file.c_str());
	if(!reader.is_open())
	{
		std::cerr << "Cannot open " << file << std::endl;
		exit(1);
	}

	std::string line;
	while(getline(reader, line))
	{
		// Read the alignments
		std::stringstream ss(line);
		std::string read1, read2;
		KAlignment ka1, ka2;	
		ss >> read1 >> ka1 >> read2 >> ka2;

		// Add the coverage of ka2 to the contig of ka1
        IDUniDataMap::iterator iter = udMap.find(ka1.contig_id);
        assert(iter != udMap.end());	
		iter->second.addOverhangCoverage(ka1.is_reverse, opt::k, ka2);
	}
    reader.close();
}

// 
// Handle command line arguments
//
void parseOptions(int argc, char** argv)
{
	bool die = false;

	// Set defaults
	opt::threshold = 50;
	opt::outfile = "unicontigs.caf";
	opt::length_cutoff = 50;
	opt::bUseDepth = 1;
	opt::bUsePairs = 0;

	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'a': arg >> opt::alignFile; break;
			case 'p': arg >> opt::pairedFile; opt::bUsePairs = 1; break;
			case 'h': arg >> opt::histFile; opt::bUsePairs = 1;  break;
			case 'l': arg >> opt::length_cutoff; break;
			case 'v': opt::verbose++; break;
			case 't': arg >> opt::threshold; break;
			case 'o': arg >> opt::outfile; break;
			case OPT_NODEPTH: opt::bUseDepth = 0; break;
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

	if(opt::threshold <= 0)
	{
		std::cerr << PROGRAM ": invalid -t,--threshold option (must be > 0)\n";
		die = true;
	}

	if(opt::length_cutoff <= 0)
	{
		std::cerr << PROGRAM ": invalid -l,--len_cutff option (must be > 0)\n";
		die = true;
	}

	if(opt::bUsePairs && opt::histFile.empty())
	{
		std::cerr << PROGRAM ": error missing --histogram file\n";
		die = true;
	}	

	if(opt::bUsePairs && opt::pairedFile.empty())
	{
		std::cerr << PROGRAM ": error missing --paired file\n";
		die = true;
	}		

	if(!opt::bUseDepth && !opt::bUsePairs)
	{
		std::cerr << PROGRAM ": invalid options, you cannot specify --no_depth without provided a --paired/hist file\n";
		die = true;
	}

	if(opt::outfile.empty())
	{
		std::cerr << PROGRAM ": invalid -o,--outfile option (cannot be empty string)\n";
		die = true;
	}

	if (argc - optind < 1) 
	{
		std::cerr << PROGRAM ": missing arguments\n";
		die = true;
	} 
	else if (argc - optind > 2) 
	{
		std::cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if(opt::alignFile.empty())
	{
		std::cerr << PROGRAM ": an alignment file must be specified\n";
		die = true;
	}

	if (die) 
	{
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
}

