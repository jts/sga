//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// OverlapThread - Threaded implementation
// of overlapping a read to a BWT/FM-index
//
#include "OverlapThread.h"

OverlapThread::OverlapThread(const OverlapAlgorithm* pOverlapper, 
							 const std::string& filename, const size_t max_items) : 
                             m_pOverlapper(pOverlapper), 
							 m_outfile(filename.c_str()), m_stopRequested(false)
{
	m_pOBList = new OverlapBlockList;
	m_pSharedWorkVec = new OverlapWorkVector;
	m_pSharedWorkVec->reserve(max_items);

	// Set up semaphores and mutexes
	sem_init( &m_producedSem, PTHREAD_PROCESS_PRIVATE, 0 );
	sem_init( &m_consumedSem, PTHREAD_PROCESS_PRIVATE, 1 );
	pthread_mutex_init(&m_mutex, NULL);
}

OverlapThread::~OverlapThread()
{
	m_outfile.close();
	delete m_pOBList;
	delete m_pSharedWorkVec;

	sem_destroy(&m_producedSem);
	sem_destroy(&m_consumedSem);
	pthread_mutex_destroy(&m_mutex);
}

void OverlapThread::swapBuffers(OverlapWorkVector* pIncoming)
{
	// Wait for the overlap thread
	// to finish with the shared vector
	// then lock the mutex
	sem_wait(&m_consumedSem);
	pthread_mutex_lock(&m_mutex);
	m_pSharedWorkVec->swap(*pIncoming);
	
	// Unlock the mutex and signal the work thread
	// that it has data
	pthread_mutex_unlock(&m_mutex);
	sem_post(&m_producedSem);
}

//
void OverlapThread::start()
{
	WARN_ONCE("check return codes");
	pthread_create(&m_thread, 0, &OverlapThread::startThread, this);
}

//
void OverlapThread::stop()
{
	// Wait for the work thread to finish
	// with the data its processing then set the stop
	// flag and unblock the work process
	sem_wait(&m_consumedSem);
	m_stopRequested = true;
	sem_post(&m_producedSem);
    pthread_join(m_thread, NULL);
}

//
void OverlapThread::run()
{
	while(1)
	{
		sem_wait(&m_producedSem);
		
		// If the stop flag was set while we waited, we are done
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
		pthread_mutex_unlock(&m_mutex);
		sem_post(&m_consumedSem);
	}
	std::cout << "Found stop signal, exitting.\n";
}

// Overlap a read
bool OverlapThread::processRead(const SeqItem& read)
{
	// Perform the work
	m_pOverlapper->overlapRead(read, m_pOBList);
	
	// Write the hits to the file
	if(!m_pOBList->empty())
	{
		WARN_ONCE("Refactor the overlap hit output code AND FIX INDICES");

		// Write the header info
		size_t numBlocks = m_pOBList->size();
		m_outfile << 0 << " " << numBlocks << " ";
		//std::cout << "<Wrote> idx: " << count << " count: " << numBlocks << "\n";
		for(OverlapBlockList::const_iterator iter = m_pOBList->begin(); iter != m_pOBList->end(); ++iter)
		{
			m_outfile << *iter << " ";
		}
		m_outfile << "\n";
		m_pOBList->clear();
	}
	return false;
}

// 
void* OverlapThread::startThread(void* obj)
{
	reinterpret_cast<OverlapThread*>(obj)->run();
	return NULL;
}
