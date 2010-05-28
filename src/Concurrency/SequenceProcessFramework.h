//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// SequenceProcessFramework - Generic framework for performing
// some operations on all sequneces in a file, serially or in parallel
//
#include "ThreadWorker.h"
#include "Timer.h"

#ifndef SEQUENCEPROCESSFRAMEWORK_H
#define SEQUENCEPROCESSFRAMEWORK_H

struct SequenceWorkItem
{
    SequenceWorkItem(size_t ri, const SeqRecord& sr) : idx(ri), read(sr) {}
    size_t idx;
    SeqRecord read;
};

namespace SequenceProcessFramework
{

// Framework for performing some processing on every sequence in a file
template<class Output, class Processor, class PostProcessor>
size_t processSequencesSerial(const std::string& readsFile, Processor* pProcessor, PostProcessor* pPostProcessor)
{
    Timer timer("SequenceProcess", true);
    SeqReader reader(readsFile);
    SeqRecord read;
    size_t currIdx = 0;
    while(reader.get(read))
    {
        SequenceWorkItem workItem(currIdx++, read);
        Output output = pProcessor->process(workItem);
        
        pPostProcessor->process(workItem, output);
        if(currIdx % 50000 == 0)
            printf("[sga] Processed %zu sequences\n", currIdx);            
    }

    double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            currIdx, proc_time_secs, (double)currIdx / proc_time_secs);    
    return currIdx;
}

// Design:
// This function is a generic function to read in sequences from
// a file and perform some work on them. The actual processing is done
// by the Processor class that is passed in. The number of threads
// created is determined by the size of the vector of processors - 
// one thread per processor. 
//
// The function reads in a batch of sequences from a file and buffers
// them. Once the buffers are full, the reads are dispatched to the thread
// which run the actual processing independently. An optional post processor
// can be specified to process the results that the threads return.
template<class Output, class Processor, class PostProcessor>
size_t processSequencesParallel(const std::string& readsFile, std::vector<Processor*> processPtrVector, PostProcessor* pPostProcessor)
{
    Timer timer("SequenceProcess", true);
    
    // The number of items to buffer before dispatching to the threads
    size_t MAX_ITEMS = 1000;
    
    // Helpful typedefs
    typedef ThreadWorker<SequenceWorkItem, Output, Processor> Thread;
    typedef std::vector<Thread*> ThreadPtrVector;

    typedef std::vector<SequenceWorkItem> SeqWorkItemVector;
    typedef std::vector<SeqWorkItemVector*> InputBufferVector;

    typedef std::vector<Output> OutputVector;
    typedef std::vector<OutputVector*> OutputBufferVector;
    typedef std::vector<sem_t*> SemaphorePtrVector;

    // Initialize threads, one thread per processor that was passed in
    int numThreads = processPtrVector.size();

    ThreadPtrVector threadVec(numThreads, NULL);
    InputBufferVector inputBuffers(numThreads, NULL);
    OutputBufferVector outputBuffers(numThreads, NULL);
    SemaphorePtrVector semVec(numThreads, NULL);

    // Create the threads
    for(int i = 0; i < numThreads; ++i)
    {
        semVec[i] = new sem_t;
        sem_init( semVec[i], PTHREAD_PROCESS_PRIVATE, 0 );

        // Create and start the thread
        threadVec[i] = new Thread(semVec[i], processPtrVector[i], MAX_ITEMS);
        threadVec[i]->start();

        inputBuffers[i] = new SeqWorkItemVector;
        inputBuffers[i]->reserve(MAX_ITEMS);

        outputBuffers[i] = new OutputVector;
        outputBuffers[i]->reserve(MAX_ITEMS);
    }

    size_t numRead = 0;
    size_t numWrote = 0;
    bool done = false;
    int next_thread = 0;
    int num_buffers_full = 0;

    SeqReader reader(readsFile);
    SeqRecord read;
    while(!done)
    {
        // Parse reads from the stream and add them into the incoming buffers
        done = !reader.get(read);
        if(!done)
        {
            inputBuffers[next_thread]->push_back(SequenceWorkItem(numRead++, read));
            if(inputBuffers[next_thread]->size() == MAX_ITEMS)
            {
                ++num_buffers_full;
                ++next_thread;
            }
        }
        
        // Once all buffers are full or the input is finished, dispatch the reads to the threads
        // by swapping work buffers. 
        if(num_buffers_full == numThreads || done)
        {
            int numLoops = 0;
            do
            {
                // Wait for all threads to be ready to receive
                for(int i = 0; i < numThreads; ++i)
                {
                    sem_wait(semVec[i]);
                    Thread* pThread = threadVec[i];
                    pThread->swapBuffers(*inputBuffers[i], *outputBuffers[i]);
                }
                num_buffers_full = 0;
                next_thread = 0;

                // Process the results
                // Write out the results of the overlaps to the asqg file
                for(int i = 0; i < numThreads; ++i)
                {
                    assert(inputBuffers[i]->size() == outputBuffers[i]->size());
                    for(size_t j = 0; j < inputBuffers[i]->size(); ++j)
                    {
                        pPostProcessor->process((*inputBuffers[i])[j], (*outputBuffers[i])[j]);
                        ++numWrote;
                    }
                    
                    inputBuffers[i]->clear();
                    outputBuffers[i]->clear();
                }

                if(numWrote % (MAX_ITEMS * numThreads * 4) == 0)
                    printf("[sga] Processed %zu sequences\n", numWrote);

                // This should never loop more than twice
                assert(numLoops < 2);
                ++numLoops;
            } while(done && numWrote < numRead);
        }
    }

    // Cleanup
    for(int i = 0; i < numThreads; ++i)
    {
        threadVec[i]->stop(); // Blocks until the thread joins
        delete threadVec[i];

        sem_destroy(semVec[i]);
        delete semVec[i];

        assert(inputBuffers[i]->empty());
        delete inputBuffers[i];

        assert(outputBuffers[i]->empty());
        delete outputBuffers[i];
    }
    
    double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            numWrote, proc_time_secs, (double)numWrote / proc_time_secs);
    return numWrote;
}

};

#endif
