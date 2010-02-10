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
	OverlapWorkItem(const SeqItem& si, bool b) : read(si), is_end(b) {}
	SeqItem read;
	bool is_end;
};


typedef std::vector<SeqItem> OverlapWorkVector;

class OverlapThread
{
	public:

		OverlapThread(const OverlapAlgorithm* pOverlapper, 
		              const std::string& filename, sem_t* pReadySem, const size_t max_items);

		~OverlapThread();

		// Exchange the contents of the shared vector with the point-to vector
		void swapBuffers(OverlapWorkVector* pIncoming);

		// Returns true if the thread is ready to receive data
		bool isReady();

		// External control functions
		void start();
		void stop();

	private:

		// Main work loop
		void run();

		// Main overlap function
		bool processRead(const SeqItem& read);
	
		// Thread entry point
		static void* startThread(void* obj);

		// Data
		
		// Private file handle and overlap list
		// which are only acccessed by the thread
		std::ofstream m_outfile;
		OverlapBlockList* m_pOBList;

		// 
		const OverlapAlgorithm* m_pOverlapper;
		pthread_t m_thread;

		// External semaphore to post to
		// when the thread is ready to receive data
		sem_t* m_pReadySem;

		// Shared data and protection variables
		sem_t m_producedSem;
		sem_t m_consumedSem;
		
		pthread_mutex_t m_mutex;
		OverlapWorkVector* m_pSharedWorkVec;
		volatile bool m_stopRequested;
		bool m_isReady;
};

#endif
