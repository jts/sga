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

#include "LockedQueue.h"
#include "Util.h"
#include "OverlapBlock.h"
#include "OverlapAlgorithm.h"

typedef LockedQueue<SeqItem> SeqItemQueue;
typedef LockedQueue<OverlapBlockList*> OBListPtrQueue;

class OverlapThread
{
	public:

		OverlapThread(const OverlapAlgorithm* pOverlapper);
		~OverlapThread();

		// 
		inline SeqItemQueue* getSeqQueue() const
		{
			return m_pSeqItemQueue;
		}

		//
		inline OBListPtrQueue* getOBListQueue() const
		{
			return m_pOBListQueue;
		}

		// External control functions
		void start();
		void stop();

	private:

		// Main work function
		void run();
	
		// Thread entry point
		static void* startThread(void* obj);

		// Data
		SeqItemQueue* m_pSeqItemQueue;
		OBListPtrQueue* m_pOBListQueue;
		const OverlapAlgorithm* m_pOverlapper;

		pthread_t m_thread;

		// The main thread will set this flag when it wants this thread to terminate
		volatile bool m_stopRequested;

};

#endif
