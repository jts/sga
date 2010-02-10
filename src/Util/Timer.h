//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Timer - Simple object to that prints the wallclock 
// length of its lifetime
//
#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <sys/time.h>

class Timer
{
	public:
		Timer(std::string s, bool silent = false) : m_desc(s), m_silent(silent)
		{
			reset();
		}

		~Timer() 
		{
			 	
			if(!m_silent) 
				printf("[Timer] %s %.2lfs\n", m_desc.c_str(), getElapsedTime()); 
		}

		double getElapsedTime() const 
		{ 
			time_t now;
			time(&now);
			return now - m_start;
		}

		void reset() { time(&m_start); }

	private:
		std::string m_desc;
		time_t m_start;
		bool m_silent;
};

#endif
