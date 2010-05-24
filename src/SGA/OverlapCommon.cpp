//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapCommon - Common wrapper used for finding overlaps
// for a set of reads
//
#include "OverlapCommon.h"

// Compute the hits for each read in the SeqReader file without threading
// Return the number of reads processed
size_t OverlapCommon::computeHitsSerial(const std::string& prefix, const std::string& readsFile, 
                                        const OverlapAlgorithm* pOverlapper, OverlapMode mode,
                                        int minOverlap, StringVector& filenameVec, std::ostream* pASQGWriter)
{
    SeqReader reader(readsFile);
    if(mode == OM_OVERLAP)
    {
        assert(pASQGWriter != NULL);
    }

    // Prepare output
    OverlapBlockList* pOutBlocks = new OverlapBlockList;

    std::string filename = prefix + HITS_EXT + GZIP_EXT;
    std::ostream* pHitsWriter = createWriter(filename);
    filenameVec.push_back(filename);

    // Overlap each read
    SeqRecord read;
    size_t currIdx = 0;
    while(reader.get(read))
    {
        if(currIdx % 50000 == 0)
            printf("[sga::overlap] Aligned %zu sequences\n", currIdx);

        if(mode == OM_OVERLAP)
        {
            OverlapResult result = pOverlapper->overlapRead(read, minOverlap, pOutBlocks);
            pOverlapper->writeResultASQG(*pASQGWriter, read, result);
            pOverlapper->writeOverlapBlocks(*pHitsWriter, currIdx++, result.isSubstring, pOutBlocks);
        }
        else
        {
            assert(mode == OM_FULLREAD);
            OverlapResult result = pOverlapper->alignReadDuplicate(read, pOutBlocks);
            pOverlapper->writeOverlapBlocks(*pHitsWriter, currIdx++, result.isSubstring, pOutBlocks);
        }
        pOutBlocks->clear();
    }

    delete pOutBlocks;
    delete pHitsWriter; 
    return currIdx;
}

// Compute the hits for each read in the SeqReader file with threading
size_t OverlapCommon::computeHitsParallel(int numThreads, const std::string& prefix,  const std::string& readsFile,  
                                          const OverlapAlgorithm* pOverlapper, OverlapMode mode,
                                          int minOverlap, StringVector& filenameVec, std::ostream* pASQGWriter)
{
    SeqReader reader(readsFile);
    if(mode == OM_OVERLAP)
    {
        assert(pASQGWriter != NULL);
    }
        
    size_t MAX_ITEMS = 1000;
    
    // Semaphore shared between all the threads indicating whether 
    // the thread is read to take data
    typedef std::vector<OverlapThread*> ThreadPtrVector;
    typedef std::vector<OverlapWorkVector*> WorkPtrVector;
    typedef std::vector<sem_t*> SemaphorePtrVector;

    // Initialize threads
    ThreadPtrVector threadVec(numThreads, NULL);
    WorkPtrVector workVec(numThreads, NULL);
    SemaphorePtrVector semVec(numThreads, NULL);

    for(int i = 0; i < numThreads; ++i)
    {
        // Create the thread
        std::stringstream ss;
        ss << prefix << "-thread" << i << HITS_EXT << GZIP_EXT;
        std::string outfile = ss.str();
        filenameVec.push_back(outfile);

        semVec[i] = new sem_t;
        sem_init( semVec[i], PTHREAD_PROCESS_PRIVATE, 0 );

        // Create and start the thread
        threadVec[i] = new OverlapThread(pOverlapper, mode, minOverlap, outfile, semVec[i], MAX_ITEMS);
        threadVec[i]->start();

        workVec[i] = new OverlapWorkVector;
        workVec[i]->reserve(MAX_ITEMS);
    }

    SeqRecord read;
    size_t numRead = 0;
    size_t numWrote = 0;
    bool done = false;
    int next_thread = 0;
    int num_buffers_full = 0;

    while(!done)
    {
        // Parse reads from the stream and add them into the incoming buffers
        done = !reader.get(read);
        if(!done)
        {
            workVec[next_thread]->push_back(OverlapWorkItem(numRead++, read));
            if(workVec[next_thread]->size() == MAX_ITEMS)
            {
                ++num_buffers_full;
                ++next_thread;
            }
        }
        
        // Once all buffers are full or the input is finished, dispatch the reads to the threads
        // by swapping work buffers. This thread then gets the OverlapResults
        // back from the worker thread which are written to the ASQG file
        // We continue in this code until all the results have been returned by the threads
        if(num_buffers_full == numThreads || done)
        {
            int numLoops = 0;
            do
            {
                // Wait for all threads to be ready to receive
                for(int i = 0; i < numThreads; ++i)
                {
                    sem_wait(semVec[i]);
                    OverlapThread* pThread = threadVec[i];
                    pThread->swapBuffers(workVec[i]);
                }
                num_buffers_full = 0;
                next_thread = 0;

                // Write out the results of the overlaps to the asqg file
                for(int i = 0; i < numThreads; ++i)
                {
                    for(size_t j = 0; j < workVec[i]->size(); ++j)
                    {
                        OverlapWorkItem& item = (*workVec[i])[j];
                        if(mode == OM_OVERLAP)
                            pOverlapper->writeResultASQG(*pASQGWriter, item.read, item.result);
                        ++numWrote;
                    }
                    workVec[i]->clear();
                }

                if(numWrote % 4000 == 0)
                    std::cout << "Aligned " << numWrote << " reads\n";

                // This should never loop more than twice
                assert(numLoops < 2);
                ++numLoops;
            } while(done && numWrote < numRead);
        }
    }

    // Stop the threads and destroy them
    for(int i = 0; i < numThreads; ++i)
    {
        threadVec[i]->stop();
        delete threadVec[i];

        sem_destroy(semVec[i]);
        delete semVec[i];

        assert(workVec[i]->empty());
        delete workVec[i];
    }

    return numWrote;
}

// Convert a line from a hits file into a vector of overlaps and sets the flag
// indicating whether the read was found to be a substring of other reads
// Only the forward read table is used since we only care about the IDs and length
// of the read, not the sequence, so that we don't need an explicit reverse read table
void OverlapCommon::parseHitsString(const std::string& hitString, 
                                    const ReadTable* pFwdRT, 
                                    const SuffixArray* pFwdSAI, const SuffixArray* pRevSAI, 
                                    size_t& readIdx, OverlapVector& outVector, bool& isSubstring)
{
    OverlapVector outvec;
    std::istringstream convertor(hitString);

    // Read the overlap blocks for a read
    size_t numBlocks;
    convertor >> readIdx >> isSubstring >> numBlocks;

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
            const ReadTable* pCurrRT = pFwdRT;
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
                if(o.id[0] < o.id[1] || (o.match.isContainment() && record.flags.isQueryRev()))
                    continue;

                outVector.push_back(o);
            }
        }
    }
}
