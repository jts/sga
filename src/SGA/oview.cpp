//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// oview - view overlap alignments
//
#include <iostream>
#include <fstream>
#include <iterator>
#include <math.h>
#include "Util.h"
#include "oview.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "LCPArray.h"
#include "SGUtil.h"
#include "OverlapQuality.h"

void detectMisalignments(const ReadTable* pRT, const OverlapMap* pOM);
void detect(const SeqItem& read, const ReadTable* pRT, const OverlapMap* pOM);


//
// Getopt
//
#define SUBPROGRAM "oview"
static const char *OVIEW_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *OVIEW_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Compute pairwise overlap between all the sequences in READS\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -i, --id=ID                      only show overlaps for read with ID\n"
"      -m, --max-overhang=D             only show D overhanging bases of the alignments (default=6)\n"
"      -c, --correct-errors             correct reads errors using overlaps\n"
"      -d, --detect-misalign            detect misaligned reads\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
	static unsigned int verbose;
	static int max_overhang;
	static std::string prefix;
	static std::string readsFile;
	static std::string readFilter;
	static bool bErrorCorrect;
	static bool bDetectMisalign;
}

static const char* shortopts = "p:m:i:cvd";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",         no_argument,       NULL, 'v' },
	{ "id",              required_argument, NULL, 'i' },
	{ "prefix",          required_argument, NULL, 'p' },
	{ "max-overhang",    required_argument, NULL, 'm' },
	{ "correct-errors",  no_argument,       NULL, 'c' },
	{ "detect-misalign", no_argument,       NULL, 'd' },
	{ "help",            no_argument,       NULL, OPT_HELP },
	{ "version",         no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int oviewMain(int argc, char** argv)
{
	parseOviewOptions(argc, argv);

	// parse reads and index them
	ReadTable* pRT = new ReadTable(opt::readsFile);
	pRT->indexReadsByID();

	// read the contain and overlap files
	std::string containFile = opt::prefix + ".ctn";
	std::string overlapFile = opt::prefix + ".ovr";
	OverlapMap overlapMap;
	parseOverlaps(containFile, overlapMap);
	parseOverlaps(overlapFile, overlapMap);


	if(opt::bDetectMisalign)
		detectMisalignments(pRT, &overlapMap);
	else if(opt::bErrorCorrect)
	{
		correctReads(pRT, &overlapMap);
	}
	else
	{
		// draw mode
		if(!opt::readFilter.empty())
		{
			drawAlignment(opt::readFilter, pRT, &overlapMap);
		}
		else
		{
			// Output each overlap
			for(size_t i = 0; i < pRT->getCount(); ++i)
				drawAlignment(pRT->getRead(i).id, pRT, &overlapMap);
		}
	}

	delete pRT;
	return 0;
}

void correctReads(const ReadTable* pRT, const OverlapMap* pOM)
{
	std::string outFile(opt::prefix + ".ec.fa");
	std::ofstream out(outFile.c_str());

	size_t total_bases = 0;
	size_t total_corrected = 0;
	for(size_t i = 0; i < pRT->getCount(); ++i)
	{
		int num_corrected = 0;

		const SeqItem& si = pRT->getRead(i);
		std::string corrected = correct(pRT->getRead(i), pRT, pOM, num_corrected);
		
		total_bases += corrected.length();
		total_corrected += num_corrected;
		out << ">" << si.id << "\n" << corrected << "\n";
	}
	printf("corrected %zu bases out of %zu total (%lf)\n", total_corrected, total_bases, (double)total_corrected / total_bases);
	out.close();
}

void detectMisalignments(const ReadTable* pRT, const OverlapMap* pOM)
{
	for(size_t i = 0; i < pRT->getCount(); ++i)
	{
		const SeqItem& read = pRT->getRead(i);
		if(read.id != "35621:20150")
			continue;

		OverlapMap::const_iterator finder = pOM->find(read.id);
		if(finder == pOM->end())
			continue;

		const OverlapVector& ov = finder->second;
		OverlapQuality qual(read, ov, pRT);

		int be_count = 0;
		for(size_t i = 0; i < ov.size(); ++i)
		{
			if(debug_getReadDistFromNames(ov[i].id[0],ov[i].id[1]) > 100)
				++be_count;
		}
		
		double improve = qual.cluster();
		printf("%d\t%lf\t%s\tBD\n", be_count, improve, read.id.c_str());
	}
}

std::string correct(const SeqItem& read, const ReadTable* pRT, const OverlapMap* pOM, int& num_corrected)
{
	OverlapMap::const_iterator finder = pOM->find(read.id);
	if(finder == pOM->end())
		return read.seq.toString();
	
	const OverlapVector & overlaps = finder->second;
	StringVector pileupVec(read.seq.length());

	for(size_t i = 0; i < read.seq.length(); ++i)
		pileupVec[i].append(1, read.seq.get(i));

	for(size_t i = 0; i < overlaps.size(); ++i)
	{
		Overlap curr = overlaps[i];

		// If the first element of the overlap is not the read
		// swap the elements
		if(curr.id[0] != read.id)
			curr.swap();

		assert(curr.id[0] == read.id);
	
		std::string otherSeq = pRT->getRead(curr.id[1]).seq.toString();

		// Make the other sequence in the same frame as the root
		if(curr.match.isRC())
		{
			otherSeq = reverseComplement(otherSeq);
			curr.match.canonize();
		}

		// Add the bases of other to the pileup
		for(size_t j = 0; j < otherSeq.length(); ++j)
		{
			// Transform the position j on otherSeq to the coordinates of the read
			int transformed = curr.match.inverseTranslate(j);
			if(transformed >= 0 && transformed < (int)pileupVec.size())
				pileupVec[transformed].append(1, otherSeq[j]);
		}
	}

	// Prior probability that there is an error at this base
	//double p_base = 0.25; // all bases are equally probable
	double p_error = 0.01;
	double threshold = 40;

	double lp_error = log(p_error);
	double lp_correct = log(1.0f - p_error);

	std::string corrected(read.seq.length(), 'A');
	num_corrected = 0;

	if(opt::verbose > 0)
		std::cout << "Pileup string for " << read.seq.toString() << "\n";

	for(size_t i = 0; i < pileupVec.size(); ++i)
	{
		double probs[4];
		std::string& str = pileupVec[i];

		// Calculate the probability of each base, assuming the sequence is haploid
		for(size_t j = 0; j < 4; ++j)
		{
			char b = ALPHABET[j];
			double lp = 0.0f;

			for(size_t k = 0; k < str.size(); ++k)
			{
				if(str[k] == b)
					lp += lp_correct;
				else
					lp += lp_error;
			}
			probs[j] = lp;
		}

		if(opt::verbose > 0)
			printf("%s A: %lf C: %lf G: %lf T: %lf\n", str.c_str(), probs[0], probs[1], probs[2], probs[3]); 

		// Correct the base if A) the most likely base is not the base in the read and B) the probability of the base is
		// above threshold
		char best = '$';
		double best_lp = 0.0f;
		for(size_t j = 0; j < 4; ++j)
		{
			if(probs[j] > best_lp || best == '$')
			{
				best = ALPHABET[j];
				best_lp = probs[j];
			}
		}

		std::sort(probs, probs+4);
		double second = probs[2];
		double diff = best_lp - second;
		(void)diff;
		(void)threshold;
		
		if(diff > threshold && best != str[0])
		{
			++num_corrected;
			corrected[i] = best;
		}
		else
		{
			corrected[i] = str[0];
		}
		//printf("best: %c bestlp: %lf second: %lf diff: %lf\n", best, best_lp, second, diff);
	}

	return corrected;
}

//
void drawAlignment(std::string rootID, const ReadTable* pRT, const OverlapMap* pOM)
{
	std::string rootSeq = pRT->getRead(rootID).seq.toString();
	DrawVector draw_vector;

	DrawData rootData(0, rootID, rootSeq, 0, 0);
	draw_vector.push_back(rootData);

	// Get all the overlaps for this read
	OverlapMap::const_iterator finder = pOM->find(rootID);
	if(finder == pOM->end())
		return;
	
	const OverlapVector& overlaps = finder->second;

	// Convert each overlap into an offset and string
	for(size_t j = 0; j < overlaps.size(); ++j)
	{
		Overlap curr = overlaps[j];

		// Swap root read into first position if necessary
		if(curr.id[0] != rootID)
			curr.swap();
		assert(curr.id[0] == rootID);

		std::string otherSeq = pRT->getRead(curr.id[1]).seq.toString();

		// 
		if(curr.match.isRC())
		{
			otherSeq = reverseComplement(otherSeq);
			curr.match.canonize();
		}
		
		// Calculate the offset between position 0 of otherSeq and the start of the root
		int offset = curr.match.inverseTranslate(0);
		draw_vector.push_back(DrawData(offset, curr.id[1], otherSeq, curr.match.getNumDiffs(), curr.match.coord[0].length()));
	}
	drawMulti(rootData.name, rootData.seq.size(), draw_vector);
}

//
void drawMulti(std::string rootName, int root_len, DrawVector& dv)
{
	if(dv.size() < 2)
		return;
	std::sort(dv.begin(), dv.end());
	int default_padding = 12;
	std::cout << "\nDrawing overlaps for read " << rootName << "\n";
	for(size_t i = 0; i < dv.size(); ++i)
	{
		int c_offset = dv[i].offset;
		int c_len = dv[i].seq.length();

		// This string runs from c_offset to c_offset + len
		// Clip the string at -max_overhang to root_len + max_overhang
		int left_clip = std::max(c_offset, -opt::max_overhang);
		int right_clip = std::min(c_offset + c_len, root_len + opt::max_overhang);
		
		// translate the clipping coordinates to the string coords
		int t_left_clip = left_clip - c_offset;
		int t_right_clip = right_clip - c_offset;
		// Calculate the length of the left padding
		int padding = default_padding + left_clip;
		std::string leader = (t_left_clip > 0) ? "..." : "";
		std::string trailer = (t_right_clip < c_len) ? "..." : ""; 
		std::string clipped = dv[i].seq.substr(t_left_clip, t_right_clip - t_left_clip);
		padding -= leader.size();
		//padding -= dv[i].name.size();
		//printf("offset: %d lc: %d rc: %d pad: %d\n", c_offset, left_clip, right_clip, padding);
		assert(padding >= 0);
		std::string padding_str(padding, ' ');
		std::string outstr = padding_str + leader + clipped + trailer;
		printf("%s\t%d\t%d\t%lf\t%s\n", outstr.c_str(), dv[i].overlapLen, dv[i].numDiff, (double)dv[i].numDiff / dv[i].overlapLen, dv[i].name.c_str());
		//std::cout /*<< dv[i].name*/ << padding_str << leader << clipped << trailer << "\t" << dv[i].name << "\n";		
	}
}

void parseOverlaps(std::string filename, OverlapMap& overlapMap)
{
	std::ifstream overlapReader(filename.c_str());
	Overlap o;
	while(overlapReader >> o)
	{
		// If we are in error correcting mode, or the read filter is not set or this read matches the filter, add it to the 
		// overlap map
		if(opt::bErrorCorrect || opt::bDetectMisalign || opt::readFilter.empty() || (o.id[0] == opt::readFilter || o.id[1] == opt::readFilter))
		{
			overlapMap[o.id[0]].push_back(o);
			overlapMap[o.id[1]].push_back(o);
		}
	}
}


// 
// Handle command line arguments
//
void parseOviewOptions(int argc, char** argv)
{
	bool die = false;

	// Defaults
	opt::max_overhang = 6;

	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'p': arg >> opt::prefix; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case 'i': arg >> opt::readFilter; break;
			case 'm': arg >> opt::max_overhang; break;
			case 'd': opt::bDetectMisalign = true; break;
			case 'c': opt::bErrorCorrect = true; break;
			case OPT_HELP:
				std::cout << OVIEW_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << OVIEW_VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 1) 
	{
		std::cerr << SUBPROGRAM ": missing arguments\n";
		die = true;
	} 
	else if (argc - optind > 1) 
	{
		std::cerr << SUBPROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) 
	{
		std::cout << "\n" << OVIEW_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::prefix.empty())
	{
		opt::prefix = stripFilename(opt::readsFile);
	}
}
