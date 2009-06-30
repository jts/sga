#include "BDE.h"
#include "Bigraph.h"
#include <math.h>
#include <iterator>

int main(int argc, char** argv)
{
	parseOptions(argc, argv);
	
	// Filenames
	std::string contigFile(argv[optind++]);
	std::string pairsFile(argv[optind++]);

	// Parse the hist and convert it to a distribution
	IntDist fragDist;
	FragmentDistribution fragSizeCounts;
	fragSizeCounts.readFromFile(opt::histFile);
	fragDist = fragSizeCounts.convertToIntDist(1.0); // Convert to a pdf and trim
	double distExpected = fragDist.expectedValue();
	std::cerr << "Fragment distribution expected value: " << distExpected << "\n";

	if(opt::verbose >= 1)
	{
		std::cerr << "Contigs file: " << contigFile << "\n";
		std::cerr << "Pairs file: " << pairsFile << "\n";
		std::cerr << "Hist file: " << opt::histFile << "\n";
		std::cerr << "Kmer: " << opt::k << "\n";
	}

	// Open distance estimate file for writing
	std::ofstream estimateWriter("project.bde");

	// Read all the contigs
	std::ifstream contigReader(contigFile.c_str());
	assert(contigReader.is_open());

	ContigMap contigMap;
	Contig c;
	while(readCAF(contigReader,c))
	{
		contigMap.insert(std::make_pair(c.getID(), c));
	}
	contigReader.close();

	PairStreamer pairStream(pairsFile);
	bool done = false;
	while(!done)
	{
		AlignPairVec apv = pairStream.getBlock();
		IDAPVecMap indexedAPVecs = splitBlock(apv);
		
		// we should always have a block of data
		if(apv.size() == 0)
		{
			done = true;
			break;
		}
	
		//
		// Process all the contigs paired to this contig
		//
		ContigID id0 = apv.front().aligns[0].contig_id;
		if(id0 != "136")
			continue;

		const Contig& contig0 = getContig(contigMap, id0);
		if(!contig0.isUnique())
			continue;
	
		ProfileVector coverageProfiles[2];
		int lowDistance = -opt::k + 1;
		int highDistance = fragDist.getEnd();
		int maxDistance = highDistance;

		for(IDAPVecMap::iterator iter = indexedAPVecs.begin(); iter != indexedAPVecs.end(); ++iter)
		{
			ContigID id1 = iter->first;
			const Contig& contig1 = getContig(contigMap, id1);

			// Skip the estimate if the contig is not unique
			if(!contig1.isUnique())
				continue;

			AlignPairVec& alignPairs = iter->second;
			std::cerr << id0 << " has " << alignPairs.size() << " alignments to " << id1 << "\n";
			std::cerr << "Lengths: " << contig0.getLength() << "," << contig1.getLength() << "\n";
			
			// Pre-process the pairs by filtering out inconsistent matches and inferring
			// the direction and relative orientation of the contigs
			EdgeDir direction;
			EdgeComp orientation;
			AlignPairVec filtered = processPairs(alignPairs, direction, orientation);

			if(direction != ED_SENSE)
				continue;
			// Calculate the provisional distances between the pairs
			// This is the sum of the distance from the start of the alignment to the 
			// end of the contig for each pair
			// Each alignment points towards the end of the contig that is in the direction
			// of the gap between the pairs
			std::vector<int> fragVec; 
			int sum = 0;
			for(AlignPairVec::const_iterator iter = filtered.begin(); iter != filtered.end(); ++iter)
			{
				int dist0 = iter->aligns[0].getDistanceToEnd(contig0.getLength());
				int dist1 = iter->aligns[1].getDistanceToEnd(contig1.getLength());
				
				//if(iter->aligns[0].align_length != iter->aligns[0].read_length || iter->aligns[1].align_length != iter->aligns[1].read_length)
				//	continue;
				//
				int prov_size = dist0 + dist1 - 2;
				fragVec.push_back(prov_size);
				sum += prov_size;
			}

			double mean = (double)sum / (double)fragVec.size();
			double expectedDistFromMean = distExpected - mean;
			std::cerr << "Simple estimate between " << id0 << " and " << id1 << " is: " << expectedDistFromMean << "\n";

			// Calculate the likelihood of each possible distance from [lowDist, highDist]
			double maxLP;
			int argmax;
			IntDist distProfile = maxPost(lowDistance, highDistance, contig0.getLength() + contig1.getLength(), fragVec, fragDist, maxLP, argmax);
			
			// Compute the coverage profile, the probability of a base at distance d being covered by the contig, from the distance profile
			int maxCoverageDist = highDistance + contig1.getLength();
			if(maxCoverageDist > maxDistance)
				maxDistance = maxCoverageDist;

			Profile covProf(id1, lowDistance, maxCoverageDist);
			for(int i = lowDistance; i <= highDistance; ++i)
			{
				//std::cout << i << "\t" << distProfile.getP(i) << "\t" << id1 << "\n";
				for(size_t j = 0; j < contig1.getLength(); ++j)
				{
					covProf.profile.addWeight(i+j, distProfile.getP(i));
				}
			}
			
			double expectation = 0;
			for(int i = lowDistance; i <= maxCoverageDist; ++i)
			{
				expectation += covProf.profile.getWeight(i);
				std::cout << i << "\t" << covProf.profile.getWeight(i) << "\t" << id1 << "\n";
			}
			//std::cerr << "Contig length " << contig1.getLength() << " expectation " << expectation << "\n";
			coverageProfiles[direction].push_back(covProf);
			//writeDistEst(estimateWriter, id0, id1, (int)expectedDistFromMean, fragVec.size());
			writeDistEst(estimateWriter, id0, id1, argmax, fragVec.size());
		}

		for(int dirIdx = 0; dirIdx < ED_COUNT; ++dirIdx)
		{
			ProfileVector& profiles = coverageProfiles[dirIdx];
			int count = 0;

			// Compute an overlap profile for each pair of coverageProfiles
			// We set the high distnace to be the longest coverage distance previously see
			// all distributions shorter than this will return 0 for the out-of-range values
			// so it is valid
			IntDist overlapProfile(lowDistance, maxDistance);
			for(ProfileVector::iterator iter = profiles.begin(); iter != profiles.end(); ++iter)
			{
				for(ProfileVector::iterator iter2 = iter + 1; iter2 != profiles.end(); ++iter2)
				{
					double expected = 0.0f;
					for(int i = lowDistance; i <= maxDistance; ++i)
					{
						double overlapP = iter->profile.getWeight(i) * iter2->profile.getWeight(i);
						expected += overlapP;
						//std::cout << i << "\t" << overlapP << "\t" << "overlap" << "\n";
						//std::cout << i << "\t" << iter->getWeight(i) << "\t" << iter2->getWeight(i) << "\t" << overlapP << "\t" << count << "\n";
					}
					std::cerr << "Expected overlap between: " << iter->id << "," << iter2->id << " = " << expected << "\n";
					++count;
				}
			}
		}
		
	}
	estimateWriter.close();
}

//
// Perform a maximum posterior estimate for the distance
//
IntDist maxPost(int minDist, int maxDist, int sumContigLens, IntVec initFrags, const IntDist& dist, double& maxLP, int& argmax)
{
	maxLP = 0;
	argmax = 0;
	bool first = true;
	
	// Set the out-of-range probability, this is represents the probability of 
	// getting a erroneous pair that is nowhere near the fragment dist
	// TODO: This is arbitrary, come up with a better way of filtering outliers
	double minP = 1.0f / 1000000.f;
	
	IntDist profile(minDist, maxDist);
	DoubleVec logPVec;

	for(int d = minDist; d <= maxDist; ++d)
	{
		// Translate the values to create new values
		// Could do the calculation in-line...
		IntVec values;
		for(IntVec::iterator iter = initFrags.begin(); iter != initFrags.end(); ++iter)
		{
			values.push_back(*iter + d);
		}

		// Make the conditional distribution
		(void)sumContigLens;
		int minFragSize = 2*opt::k + d;
		int maxFragSize = d + sumContigLens;
		IntDist conditional = makeConditional(dist, minFragSize, maxFragSize);
		//std::cout << conditional << "\n";
		int discardCount = 0;
		double lp = intVecLogP(values, conditional, minP, discardCount);
		logPVec.push_back(lp);

		//std::cout << d << "\t" << lp << "\t" << discardCount << "\n";
		if(first || lp > maxLP)
		{
			maxLP = lp;
			first = false;
			argmax = d;
		}
	}

	double sum = 0;
	for(size_t idx = 0; idx < logPVec.size(); ++idx)
	{
		sum += exp(logPVec[idx] - maxLP);
	}

	double mp = log(sum) + maxLP;
	for(size_t idx = 0; idx < logPVec.size(); ++idx)
	{
		profile.setP(minDist + idx, exp(logPVec[idx] - mp));
	}
	return profile;
}

//
// Find a contig in the map, asserting if it does not exist
//
Contig& getContig(ContigMap& cm, ContigID id)
{
	ContigMap::iterator iter = cm.find(id);
	assert(iter != cm.end());
	return iter->second;
}

//
// Write out a distance estimate
//
void writeDistEst(std::ostream& out, ContigID id0, ContigID id1, int dist, int numPairs)
{
	out << "ID0:" << id0 << "\tID1:" << id1 << "\tDE:" << dist << "\tNP:" << numPairs << std::endl;
}

//
// Calculate the log-probability of obtaining VALUES from the DIST
//
double intVecLogP(IntVec values, const IntDist& dist, double minP, int& discardCount)
{
	double sum = 0.f;
	discardCount = 0;
	for(IntVec::iterator iter = values.begin(); iter != values.end(); ++iter)
	{
		double p = dist.getP(*iter);
		if(p < minP) // avoid underflow
		{
			sum += log(minP);
			++discardCount;
		}
		else
		{
			sum += log(p);
		}
	}
	return sum;
}

//
// Make a new distribution, conditioned that the data should fall between min and max
//
IntDist makeConditional(const IntDist& original, int min, int max)
{
	IntDist out(min, max);
	for(int i = min; i <= max; ++i)
	{
		out.addWeight(i, original.getP(i));
	}
	out.normalize();
	return out;
}

//
// Filter a vector of alignment pairs to remove those that are in the
// wrong orientation
// 
AlignPairVec processPairs(const AlignPairVec& apv, EdgeDir& direction, EdgeComp& orientation)
{
	int sourceNumCanonPairs = 0;
	int sourceNumRCPairs = 0;
	AlignPairVec pairsByComp[2];
	for(AlignPairVec::const_iterator iter = apv.begin(); iter != apv.end(); ++iter)
	{
		size_t idx = iter->aligns[0].is_reverse == iter->aligns[1].is_reverse ? 1 : 0;
		if(iter->aligns[0].is_reverse)
			++sourceNumRCPairs;
		else
			++sourceNumCanonPairs;

		pairsByComp[idx].push_back(*iter);
	}

	//std::cerr << "Pairs with same comp: " << pairsByComp[1].size() << " opp comp: " << pairsByComp[0].size() << "\n";
	
	// We expect the pairs will have opposite orientation by default (Illumina pairs)
	size_t expectedDominantIdx = 1;
	size_t actualDominantIdx = (pairsByComp[0].size() > pairsByComp[1].size()) ? 0 : 1;

	// If the reads are aligned in the same direction as the contig, the pairing between
	// the contigs is in the SENSE directoin
	direction = (sourceNumCanonPairs > sourceNumRCPairs) ? ED_SENSE : ED_ANTISENSE;

	// If the read pairs are aligned between the contigs with the expected orientation, the contigs
	// have come from the same strand
	orientation = (actualDominantIdx == expectedDominantIdx) ? EC_SAME : EC_REVERSE;
	
	return pairsByComp[actualDominantIdx];
}

//
// Split a block of AlignPairs into an indexed set of AlignPairVectors
// based on the second (align[1]) id
//
IDAPVecMap splitBlock(const AlignPairVec& block)
{
	IDAPVecMap out;
	for(AlignPairVec::const_iterator iter = block.begin(); iter != block.end(); ++iter)
	{
		AlignPairVec& apv = out[iter->aligns[1].contig_id];
		apv.push_back(*iter);
	}
	return out;
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
			case 'k': arg >> opt::k; break;
			case 'h': arg >> opt::histFile; break;
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

	if (opt::k <= 0) 
	{
		std::cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}	

	if(opt::histFile.empty())
	{
		std::cerr << PROGRAM ": hist file required (use --h,--histogram)\n";
		die = true;
	}

	if (die) 
	{
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
}

