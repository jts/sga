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
#include "SequenceWorkItem.h"

#ifndef SEQUENCEPROCESSFRAMEWORK_H
#define SEQUENCEPROCESSFRAMEWORK_H

namespace SequenceProcessFramework
{

const size_t BUFFER_SIZE = 1000;

// Process n sequences from a file. With the default value of -1, n becomes the largest value representable for
// a size_t and all values will be read
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesSerial(SeqReader& reader, Processor* pProcessor, PostProcessor* pPostProcessor, size_t n = -1)
{
    Timer timer("SequenceProcess", true);

    WorkItemGenerator<Input> generator(&reader);
    Input workItem;
    
    // Generate work items using the generic generation class while the number
    // of sequences consumed from the SeqReader is less than n and there 
    // are still sequences to consume from the reader
    while(generator.getNumConsumed() < n && generator.generate(workItem))
    {
        Output output = pProcessor->process(workItem);
        
        pPostProcessor->process(workItem, output);
        if(generator.getNumConsumed() % 50000 == 0)
            printf("[sga] Processed %zu sequences\n", generator.getNumConsumed());
    }

    assert(n == (size_t)-1 || generator.getNumConsumed() == n);

    //
    double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);    
    
    return generator.getNumConsumed();
}

// Framework for performing some processing on every sequence in a file
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesSerial(const std::string& readsFile, Processor* pProcessor, PostProcessor* pPostProcessor)
{
    SeqReader reader(readsFile);
    return processSequencesSerial<Input, Output, Processor, PostProcessor>(reader, pProcessor, pPostProcessor);
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
// can be specified to process the results that the threads return. If the n
// parameter is used, at most n sequences will be read from the file
template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallel(SeqReader& reader, 
                                std::vector<Processor*> processPtrVector, 
                                PostProcessor* pPostProcessor, 
                                size_t n = -1)
{
    Timer timer("SequenceProcess", true);

    // Helpful typedefs
    typedef ThreadWorker<Input, Output, Processor> Thread;
    typedef std::vector<Thread*> ThreadPtrVector;

    typedef std::vector<Input> InputItemVector;
    typedef std::vector<InputItemVector*> InputBufferVector;

    typedef std::vector<Output> OutputVector;
    typedef std::vector<OutputVector*> OutputBufferVector;
    typedef std::vector<sem_t*> SemaphorePtrVector;

    typedef WorkItemGenerator<Input> InputGenerator;

    // Initialize threads, one thread per processor that was passed in
    int numThreads = processPtrVector.size();

    ThreadPtrVector threadVec(numThreads);
    InputBufferVector inputBuffers(numThreads);
    OutputBufferVector outputBuffers(numThreads);
    SemaphorePtrVector semVec(numThreads);

    // Create the threads
    for(int i = 0; i < numThreads; ++i)
    {
        semVec[i] = new sem_t;
        sem_init( semVec[i], PTHREAD_PROCESS_PRIVATE, 0 );

        // Create and start the thread
        threadVec[i] = new Thread(semVec[i], processPtrVector[i], BUFFER_SIZE);
        threadVec[i]->start();

        inputBuffers[i] = new InputItemVector;
        inputBuffers[i]->reserve(BUFFER_SIZE);

        outputBuffers[i] = new OutputVector;
        outputBuffers[i]->reserve(BUFFER_SIZE);
    }

    size_t numWorkItemsRead = 0;
    size_t numWorkItemsWrote = 0;
    bool done = false;
    int next_thread = 0;
    int num_buffers_full = 0;

    InputGenerator generator(&reader);
    while(!done)
    {
        // Parse reads from the stream and add them into the incoming buffers
        Input workItem;
        bool valid = generator.generate(workItem);
        if(valid)
        {
            inputBuffers[next_thread]->push_back(workItem);
            numWorkItemsRead += 1;

            // Change buffers if this one is full
            if(inputBuffers[next_thread]->size() == BUFFER_SIZE)
            {
                ++num_buffers_full;
                ++next_thread;
            }
        }
        
        done = !valid || generator.getNumConsumed() == n;

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

                // Process the results and clear the buffers
                for(int i = 0; i < numThreads; ++i)
                {
                    assert(inputBuffers[i]->size() == outputBuffers[i]->size());
                    for(size_t j = 0; j < inputBuffers[i]->size(); ++j)
                    {
                        pPostProcessor->process((*inputBuffers[i])[j], (*outputBuffers[i])[j]);
                        ++numWorkItemsWrote;
                    }
                    
                    inputBuffers[i]->clear();
                    outputBuffers[i]->clear();
                }

                if(generator.getNumConsumed() % (50 * BUFFER_SIZE * numThreads) == 0)
                    printf("[sga] Processed %zu sequences\n", generator.getNumConsumed());

                // This should never loop more than twice
                assert(numLoops < 2);
                ++numLoops;
            } while(done && numWorkItemsWrote < numWorkItemsRead);
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
    assert(n == (size_t)-1 || generator.getNumConsumed() == n);
    assert(numWorkItemsRead == numWorkItemsWrote);

    double proc_time_secs = timer.getElapsedWallTime();
    printf("[sga::process] processed %zu sequences in %lfs (%lf sequences/s)\n", 
            generator.getNumConsumed(), proc_time_secs, (double)generator.getNumConsumed() / proc_time_secs);
    return generator.getNumConsumed();
}

template<class Input, class Output, class Processor, class PostProcessor>
size_t processSequencesParallel(const std::string& readsFile, std::vector<Processor*> processPtrVector, PostProcessor* pPostProcessor)
{
    SeqReader reader(readsFile);
    return processSequencesParallel<Input, Output, Processor, PostProcessor>(reader, processPtrVector, pPostProcessor);
}

};

#endif
