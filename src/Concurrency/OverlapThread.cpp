//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// OverlapThread - Threaded implementation
// of overlapping a read to a BWT/FM-index
//
#include "OverlapThread.h"

//
void OverlapThread::start()
{
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
		std::cout << "OverlapThread::work()\n";
		sleep(1);
	}
}

// 
void* OverlapThread::startThread(void* obj)
{
	reinterpret_cast<OverlapThread*>(obj)->run();
	return NULL;
}
