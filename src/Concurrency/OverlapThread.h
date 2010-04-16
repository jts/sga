//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// OverlapThread - Threaded function
// for aligning a read to a BWT/FM-Index
//
#ifndef OVERLAPTHREAD_H
#define OVERLAPTHREAD_H

#include <semaphore.h>
#include "Util.h"
#include "OverlapBlock.h"
#include "OverlapAlgorithm.h"

struct OverlapWorkItem
{
	OverlapWorkItem(size_t ri, const SeqRecord& sr) : idx(ri), read(sr) {}
	size_t idx;
	SeqRecord read;
	OverlapResult result;
};


typedef std::vector<OverlapWorkItem> OverlapWorkVector;

class OverlapThread
{
	public:

		OverlapThread(const OverlapAlgorithm* pOverlapper, 
		              const std::string& filename, sem_t* pReadySem, const size_t max_items);

		~OverlapThread();

		// Exchange the contents of the shared vector with the point-to vector
		void swapBuffers(OverlapWorkVector* pIncoming);

		// External control functions
		void start();
		void stop();
		bool isReady();

	private:

		// Main work loop
		void run();

		// Main overlap function
		void processRead(OverlapWorkItem& read);
	
		// Thread entry point
		static void* startThread(void* obj);

		
		// Private file handle and overlap list
		// which are only acccessed by the thread
		std::ostream* m_pOutfile;
		OverlapBlockList* m_pOBList;

		// Handles
		const OverlapAlgorithm* m_pOverlapper;
		pthread_t m_thread;

		// External semaphore to post to
		// when the thread is ready to receive data
		sem_t* m_pReadySem;

		// Semaphore the external caller posts to when
		// data is ready to consume
		sem_t m_producedSem;
		
		// Shared data
		pthread_mutex_t m_mutex;
		OverlapWorkVector* m_pSharedWorkVec;
		volatile bool m_stopRequested;
		bool m_isReady;
};

#endif
