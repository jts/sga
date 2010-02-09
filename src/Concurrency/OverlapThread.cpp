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

OverlapThread::OverlapThread(const OverlapAlgorithm* pOverlapper) : m_pOverlapper(pOverlapper), m_stopRequested(false)
{
	m_pSeqItemQueue = new SeqItemQueue;
	m_pOBListQueue = new OBListPtrQueue;
}

OverlapThread::~OverlapThread()
{
	delete m_pSeqItemQueue;
	delete m_pOBListQueue;
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
	// Signal thread to stop and wait for it to rejoin
	m_stopRequested = true;
    pthread_join(m_thread, NULL);
}

//
void OverlapThread::run()
{
	while(!m_stopRequested)
	{
		if(!m_pSeqItemQueue->empty())
		{
			// Take an element from the queue
			SeqItem read = m_pSeqItemQueue->pop();

			// We allocate a list here to store the result of the overlap
			// which is deallocated by the main thread
			OverlapBlockList* pList = new OverlapBlockList;
			m_pOverlapper->overlapRead(read, pList);
			m_pOBListQueue->push(pList);
		}
	}

	std::cout << "Thread exiting\n";
}

// 
void* OverlapThread::startThread(void* obj)
{
	reinterpret_cast<OverlapThread*>(obj)->run();
	return NULL;
}
