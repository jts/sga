//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// LockedQueue.h - Wrapper for the STL
// queue that uses mutexes to ensure
// concurrency of reads/writes
//
#ifndef LOCKEDQUEUE_H
#define LOCKEDQUEUE_H

#include <queue>
#include <pthread.h>
#include <semaphore.h>
#include "Util.h"

template<class T>
class LockedQueue
{
	public:
	
		LockedQueue()
		{
			pthread_mutex_init(&m_mutex, NULL);
			size_t MAX_ITEMS = 10;
			sem_init( &m_emptySem, PTHREAD_PROCESS_PRIVATE, 0 );
			sem_init( &m_fullSem, PTHREAD_PROCESS_PRIVATE, MAX_ITEMS );

		}

		~LockedQueue()
		{
			pthread_mutex_destroy(&m_mutex);
			sem_destroy(&m_emptySem);
			sem_destroy(&m_fullSem);
		}

		// Pop differs here from the STL queue as it returns (a copy) of the first element
		T pop()
		{
			// Wait on the empty-ness semaphore
			// then pull a single item from the queue
			// and post the full-ness semaphore
			sem_wait(&m_emptySem);
			pthread_mutex_lock(&m_mutex);
			T data = m_queue.front();
			m_queue.pop();
			sem_post(&m_fullSem);
			pthread_mutex_unlock(&m_mutex);
			return data;
		}

		void push(const T& data)
		{
			// Wait until the queue is not full then 
			// lock the queue, add the data and increase 
			// the count of the empty-ness semaphore
			sem_wait(&m_fullSem);
			pthread_mutex_lock(&m_mutex);
			m_queue.push(data);
			sem_post(&m_emptySem);
			pthread_mutex_unlock(&m_mutex);
		}

		bool empty()
		{
			pthread_mutex_lock(&m_mutex);
			bool b = m_queue.empty();
			pthread_mutex_unlock(&m_mutex);
			return b;
		}

	private:

		std::queue<T> m_queue;
		sem_t m_emptySem;
		sem_t m_fullSem;
		pthread_mutex_t m_mutex;
};

#endif
