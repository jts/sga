//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// OverlapThread - Threaded implementation
// of overlapping a read to a BWT/FM-index
//
#ifndef OVERLAPTHREAD_H
#define OVERLAPTHREAD_H

typedef LockedQueue<SeqItem> SeqItemQueue;
typedef LockedQueue<OverlapBlockList*> OBListLQueue;

class OverlapThread
{
	public:

		OverlapThread(SeqItemQueue* pSIQ, OBListQueue* pOBLQ) :
		              m_pSeqItemQueue(pSIQ), m_pOBListQueue(pOBLQ) {}

		void start();
		void run();
		void stop();

	private:

		SeqItemQueue m_pSeqItemQueue;
		OBListQueue m_pOBListQueue;
};

#endif
