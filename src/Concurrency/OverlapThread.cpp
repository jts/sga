///-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapThread - Threaded implementation
// of overlapping a read to a BWT/FM-index
//
#include "OverlapThread.h"

//
OverlapThread::OverlapThread(const OverlapAlgorithm* pOverlapper, 
							 const std::string& filename, 
							 sem_t* pReadySem, 
							 const size_t max_items) : m_outfile(filename.c_str()), 
							 m_pOverlapper(pOverlapper),
							 m_pReadySem(pReadySem), 
							 m_stopRequested(false), 
							 m_isReady(false)
{
	m_pOBList = new OverlapBlockList;
	m_pSharedWorkVec = new OverlapWorkVector;
	m_pSharedWorkVec->reserve(max_items);

	// Set up semaphores and mutexes
	sem_init( &m_producedSem, PTHREAD_PROCESS_PRIVATE, 0 );
	pthread_mutex_init(&m_mutex, NULL);
}

//
OverlapThread::~OverlapThread()
{
	m_outfile.close();
	delete m_pOBList;
	delete m_pSharedWorkVec;

	sem_destroy(&m_producedSem);
	pthread_mutex_destroy(&m_mutex);
}

// Externally-called function from the main thread
// It locks the mutex and swaps the vectors, providing
// data for the work thread to operate on
void OverlapThread::swapBuffers(OverlapWorkVector* pIncoming)
{
	pthread_mutex_lock(&m_mutex);
	
	m_isReady = false;
	m_pSharedWorkVec->swap(*pIncoming);
	
	pthread_mutex_unlock(&m_mutex);
	sem_post(&m_producedSem);
}

// Externally-called function to start the worker
void OverlapThread::start()
{
	WARN_ONCE("check return codes");
	pthread_create(&m_thread, 0, &OverlapThread::startThread, this);
}

// Externally-called function to tell the worker to stop
void OverlapThread::stop()
{
	// Wait for the work thread to finish
	// with the data its processing then set the stop
	// flag and unblock the work process
	sem_wait(m_pReadySem);
	m_stopRequested = true;
	sem_post(&m_producedSem);
    pthread_join(m_thread, NULL);
}

// 
bool OverlapThread::isReady()
{
	return m_isReady;
}

// Main worker loop
void OverlapThread::run()
{
	// Indicate that the thread is ready to receive data
	pthread_mutex_lock(&m_mutex);
	m_isReady = true;
	sem_post(m_pReadySem);
	pthread_mutex_unlock(&m_mutex);

	while(1)
	{
		// Block until there is some data to use
		sem_wait(&m_producedSem);
		
		// If the stop flag is now set, finish
		if(m_stopRequested)
			break;

		// Lock the shared buffer and process all
		// the contained reads
		pthread_mutex_lock(&m_mutex);

		size_t num_reads = m_pSharedWorkVec->size();
		for(size_t i = 0; i < num_reads; ++i)
		{
			processRead((*m_pSharedWorkVec)[i]);
		}
		m_pSharedWorkVec->clear();
		
		m_isReady = true;
		pthread_mutex_unlock(&m_mutex);

		// Post to the semaphores
		sem_post(m_pReadySem);
	}
}

// Overlap a read
bool OverlapThread::processRead(const OverlapWorkItem& item)
{
	m_pOverlapper->overlapRead(item.read, m_pOBList);
	m_pOverlapper->writeOverlapBlocks(item.idx, m_pOBList, m_outfile);
	m_pOBList->clear();
	return false;
}

// Thread entry point
void* OverlapThread::startThread(void* obj)
{
	reinterpret_cast<OverlapThread*>(obj)->run();
	return NULL;
}
