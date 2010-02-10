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
#include "Util.h"

template<class T>
class LockedQueue
{
	public:
	
		LockedQueue()
		{
			pthread_mutex_init(&m_mutex, NULL);
		}

		~LockedQueue()
		{
			pthread_mutex_destroy(&m_mutex);
		}

		// Pop differs here from the STL queue as it returns (a copy) of the first element
		T pop()
		{
			pthread_mutex_lock(&m_mutex);
			T data = m_queue.front();
			m_queue.pop();
			pthread_mutex_unlock(&m_mutex);
			return data;
		}

		void push(const T& data)
		{
			pthread_mutex_lock(&m_mutex);
			m_queue.push(data);
			pthread_mutex_unlock(&m_mutex);
		}

		size_t size()
		{
			pthread_mutex_lock(&m_mutex);
			size_t n = m_queue.size();
			pthread_mutex_unlock(&m_mutex);
			return n;
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
		pthread_mutex_t m_mutex;
};

#endif
