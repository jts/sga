//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// overlap - compute pairwise overlaps between reads
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "overlap.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "LCPArray.h"
#include "SGACommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "AssembleExact.h"
#include "OverlapThread.h"

//
// Getopt
//
#define SUBPROGRAM "overlap"
static const char *OVERLAP_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *OVERLAP_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Compute pairwise overlap between all the sequences in READS\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads to compute the overlaps (default: 1)\n"
"      -e, --error-rate                 the maximum error rate allowed to consider two sequences aligned\n"
"      -m, --min-overlap=LEN            minimum overlap required between two reads\n"
"      -p, --prefix=PREFIX              use PREFIX instead of the prefix of the reads filename for the input/output files\n"
"      -d, --max-diff=D                 report all prefix-suffix matches that have at most D differences\n"
"      -i, --irreducible                only output the irreducible edges for each node\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
	static unsigned int verbose;
	static unsigned int minOverlap = DEFAULT_MIN_OVERLAP;
	static unsigned int maxDiff = 0;
	static int numThreads = 1;
	static double errorRate;
	static std::string prefix;
	static std::string readsFile;
	static bool bIrreducibleOnly;
}

static const char* shortopts = "p:m:d:e:t:vi";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "threads",     required_argument, NULL, 't' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "max-diff",    required_argument, NULL, 'd' },
	{ "prefix",      required_argument, NULL, 'p' },
	{ "error-rate",  required_argument, NULL, 'e' },
	{ "irreducible", no_argument,       NULL, 'i' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int overlapMain(int argc, char** argv)
{
	Timer* pTimer = new Timer(PROGRAM_IDENT);
	parseOverlapOptions(argc, argv);
	
	StringVector hitsFilenames;
	computeHitsBWT(hitsFilenames);
	convertHitsToOverlaps(hitsFilenames);

	delete pTimer;
	if(opt::numThreads > 1)
		pthread_exit(NULL);

	return 0;
}

void computeHitsBWT(StringVector& filenameVec)
{
	BWT* pBWT = new BWT(opt::prefix + BWT_EXT);
	BWT* pRBWT = new BWT(opt::prefix + RBWT_EXT);
	OverlapAlgorithm* pOverlapper = new OverlapAlgorithm(pBWT, pRBWT, opt::minOverlap, opt::errorRate, opt::bIrreducibleOnly);
	SeqReader reader(opt::readsFile);

	// Start timing
	Timer timer("BWT Alignment", true);

	size_t count = 0;
	if(opt::numThreads <= 1)
	{
		count = computeHitsSerial(reader, pOverlapper, filenameVec);
	}
	else
	{
		count = computeHitsParallelBatch(reader, pOverlapper, filenameVec);
	}

	double align_time_secs = timer.getElapsedWallTime();
	printf("[%s] aligned %zu sequences in %lfs (%lf sequences/s)\n", PROGRAM_IDENT, count, align_time_secs, (double)count / align_time_secs);

	// 
	delete pBWT;
	delete pRBWT;
	delete pOverlapper;
}

// Compute the hits for each read in the SeqReader file without threading
// Return the number of reads processed
size_t computeHitsSerial(SeqReader& reader, const OverlapAlgorithm* pOverlapper, StringVector& filenameVec)
{
	printf("[%s] starting serial-mode overlap computation\n", PROGRAM_IDENT);

	OverlapBlockList* pOutBlocks = new OverlapBlockList;

	std::string filename = opt::prefix + HITS_EXT;
	std::ofstream writer(filename.c_str());
	assert(writer.is_open());
	filenameVec.push_back(filename);

	SeqRecord read;
	size_t currIdx = 0;
	while(reader.get(read))
	{
		if(opt::verbose > 0 && currIdx % 50000 == 0)
			printf("[%s] Aligned %zu sequences\n", PROGRAM_IDENT, currIdx);

		pOverlapper->overlapRead(read, pOutBlocks);
		pOverlapper->writeOverlapBlocks(currIdx++, pOutBlocks, writer);
		pOutBlocks->clear();
	}

	delete pOutBlocks;
	writer.close();

	return currIdx;
}

// Compute the hits for each read in the SeqReader file with threading
// Design: The worker threads are encapsulated in the OverlapThread class
// The main thread reads in data and dispatches to a waiting worker thread.
// All threads share a semaphore indicating whether they can accept data or not
// The main thread and each worker share another semaphore which blocks the worker
// until data has been dispatched to the worker by main.
// The total number of reads processed is returned
size_t computeHitsParallel(SeqReader& reader, const OverlapAlgorithm* pOverlapper, StringVector& filenameVec)
{
	printf("[%s] starting parallel-mode overlap computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);

	size_t MAX_ITEMS = 100;
	
	// Semaphore shared between all the threads indicating whether 
	// the thread is read to take data
	sem_t readySem;
	sem_init( &readySem, PTHREAD_PROCESS_PRIVATE, 0 );

	typedef std::vector<OverlapThread*> ThreadPtrVector;
	typedef std::queue<OverlapWorkVector*> WorkPtrList;
	
	OverlapWorkVector* pWorkVector = new OverlapWorkVector;
	pWorkVector->reserve(MAX_ITEMS);
	
	WorkPtrList freeQueue;
	WorkPtrList fullQueue;

	// Initialize threads
	ThreadPtrVector threadVec(opt::numThreads, NULL);
	for(int i = 0; i < opt::numThreads; ++i)
	{
		// Create the thread
		std::stringstream ss;
		ss << opt::prefix << "-thread" << i << HITS_EXT;
		std::string outfile = ss.str();
		filenameVec.push_back(outfile);

		// Create and start the thread
		threadVec[i] = new OverlapThread(pOverlapper, outfile, &readySem, MAX_ITEMS);
		threadVec[i]->start();

		OverlapWorkVector* pWorkVec = new OverlapWorkVector;
		pWorkVec->reserve(MAX_ITEMS);
		freeQueue.push(pWorkVec);
	}

	SeqRecord read;
	size_t currIdx = 0;
	bool done = false;
	
	while(!done)
	{
		assert(!freeQueue.empty());
		// Parse one read from the stream
		done = !reader.get(read);
		if(!done)
		{
			freeQueue.front()->push_back(OverlapWorkItem(currIdx++, read));
			if(freeQueue.front()->size() == MAX_ITEMS)
			{
				fullQueue.push(freeQueue.front());
				freeQueue.pop();
			}
		}
		else
		{
			// Move all the buffers to the full list
			while(!freeQueue.empty())
			{
				fullQueue.push(freeQueue.front());
				freeQueue.pop();
			}
		}

		// All the buffers are full, dispatch the work to the threads
		if(!fullQueue.empty())
		{			
			// Block on the semaphore if there are no free queues
			// to read data into, otherwise just do trywait
			int code = 0;
			while(!fullQueue.empty() && (code == 0 || done))
			{
				if(freeQueue.empty())
					code = sem_wait(&readySem);
				else
					code = sem_trywait(&readySem);
				
				// sem_wait/trywait returns 0 if the semaphore was successfully 
				// acquired and decremented
				if(code == 0)
				{
					// Determine which thread is now available and dispatch the
					// data to it
					int selected_id = -1;
					for(int i = 0; i < opt::numThreads; ++i)
					{
						if(threadVec[i]->isReady())
						{
							selected_id = i;
							break;
						}
					}

					assert(selected_id >= 0);
					OverlapThread* pThread = threadVec[selected_id];
					OverlapWorkVector* pCurrBuffer = fullQueue.front();
					fullQueue.pop();

					pThread->swapBuffers(pCurrBuffer);
					assert(pCurrBuffer->empty());
					freeQueue.push(pCurrBuffer);
				}
			}
		}
		
		/*
		// All the buffers are full, dispatch the work to the threads
		if(freeQueue.empty() || done)
		{
			while(!fullQueue.empty())
			{
				// Wait for one of the work threads to become available
				sem_wait(&readySem);

				// Determine which thread is now available and dispatch the
				// data to it
				int selected_id = -1;
				for(int i = 0; i < opt::numThreads; ++i)
				{
					if(threadVec[i]->isReady())
					{
						selected_id = i;
						break;
					}
				}

				assert(selected_id >= 0);
				OverlapThread* pThread = threadVec[selected_id];
				OverlapWorkVector* pCurrBuffer = fullQueue.front();
				fullQueue.pop();


				pThread->swapBuffers(pCurrBuffer);
				assert(pCurrBuffer->empty());
				freeQueue.push(pCurrBuffer);
			}
		}
		*/
	}
	assert(pWorkVector->empty());
	

	/*
	while(!done)
	{
		// Parse one read from the stream
		done = !reader.get(read);
		if(!done)
			pWorkVector->push_back(OverlapWorkItem(currIdx++, read));

		// The buffer is full or there are no more 
		// sequences to read, push the contents of the 
		// buffer to a work thread
		if(pWorkVector->size() == MAX_ITEMS || done)
		{
			// Wait for one of the work threads to become available
			sem_wait(&readySem);

			// Determine which thread is now available and dispatch the
			// data to it
			int selected_id = -1;
			for(int i = 0; i < opt::numThreads; ++i)
			{
				if(threadVec[i]->isReady())
				{
					selected_id = i;
					break;
				}
			}
			assert(selected_id >= 0);

			OverlapThread* pThread = threadVec[selected_id];
			
			// Dispatch the data to the work thread which clears the buffer
			pThread->swapBuffers(pWorkVector);
			assert(pWorkVector->empty());
		}
	}
	assert(pWorkVector->empty());
	*/

	// Stop the threads and destroy them
	for(int i = 0; i < opt::numThreads; ++i)
	{
		threadVec[i]->stop();
		delete threadVec[i];
	}

	assert(fullQueue.empty());
	while(!freeQueue.empty())
	{
		OverlapWorkVector* pCurr = freeQueue.front();
		freeQueue.pop();
		delete pCurr;
	}
	sem_destroy(&readySem);
	delete pWorkVector;
	return currIdx;
}

// Compute the hits for each read in the SeqReader file with threading
size_t computeHitsParallelBatch(SeqReader& reader, const OverlapAlgorithm* pOverlapper, StringVector& filenameVec)
{
	printf("[%s] starting parallel-mode overlap computation with %d threads\n", PROGRAM_IDENT, opt::numThreads);
	size_t MAX_ITEMS = 1000;
	
	// Semaphore shared between all the threads indicating whether 
	// the thread is read to take data
	typedef std::vector<OverlapThread*> ThreadPtrVector;
	typedef std::vector<OverlapWorkVector*> WorkPtrVector;
	typedef std::vector<sem_t*> SemaphorePtrVector;

	// Initialize threads
	ThreadPtrVector threadVec(opt::numThreads, NULL);
	WorkPtrVector workVec(opt::numThreads, NULL);
	SemaphorePtrVector semVec(opt::numThreads, NULL);

	for(int i = 0; i < opt::numThreads; ++i)
	{
		// Create the thread
		std::stringstream ss;
		ss << opt::prefix << "-thread" << i << HITS_EXT;
		std::string outfile = ss.str();
		filenameVec.push_back(outfile);

		semVec[i] = new sem_t;
		sem_init( semVec[i], PTHREAD_PROCESS_PRIVATE, 0 );

		// Create and start the thread
		threadVec[i] = new OverlapThread(pOverlapper, outfile, semVec[i], MAX_ITEMS);
		threadVec[i]->start();

		workVec[i] = new OverlapWorkVector;
		workVec[i]->reserve(MAX_ITEMS);
	}

	SeqRecord read;
	size_t currIdx = 0;
	bool done = false;
	int next_thread = 0;
	int	num_buffers_full = 0;

	while(!done)
	{
		// Parse one read from the stream
		done = !reader.get(read);
		if(!done)
		{
			workVec[next_thread]->push_back(OverlapWorkItem(currIdx++, read));
			if(workVec[next_thread]->size() == MAX_ITEMS)
			{
				++num_buffers_full;
				++next_thread;
			}
		}

		if(num_buffers_full == opt::numThreads || done)
		{
			// Wait for all threads to be ready to receive
			for(int i = 0; i < opt::numThreads; ++i)
			{
				sem_wait(semVec[i]);
				OverlapThread* pThread = threadVec[i];
				pThread->swapBuffers(workVec[i]);
				assert(workVec[i]->empty());
			}
			num_buffers_full = 0;
			next_thread = 0;
		}
	}

	// Stop the threads and destroy them
	for(int i = 0; i < opt::numThreads; ++i)
	{
		threadVec[i]->stop();
		delete threadVec[i];

		sem_destroy(semVec[i]);
		delete semVec[i];

		assert(workVec[i]->empty());
		delete workVec[i];
	}

	return currIdx;
}



//
void convertHitsToOverlaps(const StringVector& hitsFilenames)
{
	printf("[%s] converting suffix array interval hits to overlaps\n", PROGRAM_IDENT);

	// Load the suffix array index and the reverse suffix array index
	// Note these are not the full suffix arrays
	SuffixArray* pFwdSAI = new SuffixArray(opt::prefix + SAI_EXT);
	SuffixArray* pRevSAI = new SuffixArray(opt::prefix + RSAI_EXT);

	// Load the read tables
	ReadTable* pFwdRT = new ReadTable(opt::readsFile);
	ReadTable* pRevRT = new ReadTable();
	pRevRT->initializeReverse(pFwdRT);

	// Open files output files
	std::string overlapFile = opt::prefix + OVR_EXT;
	std::ofstream overlapHandle(overlapFile.c_str());
	assert(overlapHandle.is_open());

	std::string containFile = opt::prefix + CTN_EXT;
	std::ofstream containHandle(containFile.c_str());
	assert(containHandle.is_open());

	for(StringVector::const_iterator iter = hitsFilenames.begin(); iter != hitsFilenames.end(); ++iter)
	{
		printf("[%s] parsing file %s\n", PROGRAM_IDENT, iter->c_str());
		std::ifstream reader(iter->c_str());
	
		// Read each hit sequentially, converting it to an overlap
		std::string line;
		while(getline(reader, line))
		{
			std::istringstream convertor(line);

			// Read the overlap blocks for a read
			size_t readIdx;
			size_t numBlocks;
			convertor >> readIdx >> numBlocks;

			//std::cout << "<Read> idx: " << readIdx << " count: " << numBlocks << "\n";
			for(size_t i = 0; i < numBlocks; ++i)
			{
				// Read the block
				OverlapBlock record;
				convertor >> record;
				//std::cout << "\t" << record << "\n";

				// Iterate through the range and write the overlaps
				for(int64_t j = record.ranges.interval[0].lower; j <= record.ranges.interval[0].upper; ++j)
				{
					const ReadTable* pCurrRT = (record.flags.isTargetRev()) ? pRevRT : pFwdRT;
					const SuffixArray* pCurrSAI = (record.flags.isTargetRev()) ? pRevSAI : pFwdSAI;
					const SeqItem& query = pCurrRT->getRead(readIdx);

					int64_t saIdx = j;

					// The index of the second read is given as the position in the SuffixArray index
					const SeqItem& target = pCurrRT->getRead(pCurrSAI->get(saIdx).getID());

					// Skip self alignments and non-canonical (where the query read has a lexo. higher name)
					if(query.id != target.id)
					{	
						// Compute the endpoints of the overlap
						int s1 = query.seq.length() - record.overlapLen;
						int e1 = s1 + record.overlapLen - 1;
						SeqCoord sc1(s1, e1, query.seq.length());

						int s2 = 0; // The start of the second hit must be zero by definition of a prefix/suffix match
						int e2 = s2 + record.overlapLen - 1;
						SeqCoord sc2(s2, e2, target.seq.length());

						// The coordinates are always with respect to the read, so flip them if
						// we aligned to/from the reverse of the read
						if(record.flags.isQueryRev())
							sc1.flip();
						if(record.flags.isTargetRev())
							sc2.flip();

						bool isRC = record.flags.isTargetRev() != record.flags.isQueryRev();

						Overlap o(query.id, sc1, target.id, sc2, isRC, record.numDiff);
					
						// The alignment logic above has the potential to produce duplicate alignments
						// To avoid this, we skip overlaps where the id of the first coord is lexo. lower than 
						// the second or the match is a containment and the query is reversed (containments can be 
						// output up to 4 times total).
						// If we are running in irreducible mode this is not necessary
						if(o.id[0] < o.id[1] || (o.match.isContainment() && record.flags.isQueryRev()))
							continue;
						writeOverlap(o, containHandle, overlapHandle);
					}
				}
			}
		}
	}

	// Delete allocated data
	delete pFwdSAI;
	delete pRevSAI;
	delete pFwdRT;
	delete pRevRT;

	// Close files
	overlapHandle.close();
	containHandle.close();
}

// Before sanity checks on the overlaps and write them out
void writeOverlap(Overlap& ovr, std::ofstream& containHandle, std::ofstream& overlapHandle)
{
	// Ensure that the overlap is not a containment
	if(ovr.match.coord[0].isContained() || ovr.match.coord[1].isContained())
	{
		containHandle << ovr << "\n";
		return;
	}

	// Unless both overlaps are extreme, skip
	if(!ovr.match.coord[0].isExtreme() || !ovr.match.coord[1].isExtreme())
	{
		std::cerr << "Skipping non-extreme overlap: " << ovr << "\n";
		return;
	}

	bool sameStrand = !ovr.match.isRC();
	bool proper = false;
	
	if(sameStrand)
	{
		proper = (ovr.match.coord[0].isLeftExtreme() != ovr.match.coord[1].isLeftExtreme() && 
				  ovr.match.coord[0].isRightExtreme() != ovr.match.coord[1].isRightExtreme());
	}
	else
	{
		proper = (ovr.match.coord[0].isLeftExtreme() == ovr.match.coord[1].isLeftExtreme() && 
				  ovr.match.coord[0].isRightExtreme() == ovr.match.coord[1].isRightExtreme());
	}
	
	if(!proper)
	{
		std::cerr << "Skipping improper overlap: " << ovr << "\n";
		return;
	}

	// All checks passed, output the overlap
	overlapHandle << ovr << "\n";
}

// 
// Handle command line arguments
//
void parseOverlapOptions(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) 
		{
			case 'm': arg >> opt::minOverlap; break;
			case 'p': arg >> opt::prefix; break;
			case 'e': arg >> opt::errorRate; break;
			case 'd': arg >> opt::maxDiff; break;
			case 't': arg >> opt::numThreads; break;
			case 'i': opt::bIrreducibleOnly = true; break;
			case '?': die = true; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cout << OVERLAP_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cout << OVERLAP_VERSION_MESSAGE;
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

	if(opt::numThreads <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
		die = true;
	}

	if (die) 
	{
		std::cout << "\n" << OVERLAP_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}

	// Parse the input filenames
	opt::readsFile = argv[optind++];

	if(opt::prefix.empty())
	{
		opt::prefix = stripFilename(opt::readsFile);
	}
}
