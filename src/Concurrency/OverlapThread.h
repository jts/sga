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
#include "LockedQueue.h"
#include "Util.h"
#include "OverlapBlock.h"
#include "OverlapAlgorithm.h"

struct OverlapWorkItem
{
	OverlapWorkItem(const SeqItem& si, bool b) : read(si), is_end(b) {}
	SeqItem read;
	bool is_end;
};

typedef LockedQueue<OverlapWorkItem> OverlapWorkQueue;

class OverlapThread
{
	public:

		OverlapThread(const OverlapAlgorithm* pOverlapper, const std::string& filename);
		~OverlapThread();

		// 
		void add(const SeqItem& read);

		// 
		inline OverlapWorkQueue* getWorkQueue() const
		{
			return m_pWorkQueue;
		}

		void waitConsumed();
		void postProduced();
		void lockMutex();
		void unlockMutex();

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
		const OverlapAlgorithm* m_pOverlapper;

		// Incoming queue, protected by a mutex
		OverlapWorkQueue* m_pWorkQueue;

		// Private file handle and overlap list
		// which are only acccessed by the thread
		std::ofstream m_outfile;
		OverlapBlockList* m_pOBList;

		// 
		pthread_t m_thread;

		// Shared data and protection variables
		sem_t m_producedSem;
		sem_t m_consumedSem;
		pthread_mutex_t m_mutex;
		SeqItem m_sharedSeq;
		std::vector<SeqItem> m_sharedSeqVec;

		static const int MAX_NUM_ITEMS = 50000;

		// The main thread will set this flag when it wants this thread to terminate
		volatile bool m_stopRequested;

};

#endif
