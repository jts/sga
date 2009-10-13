//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Timer - Simple object to that prints the length of its lifetime
//
#ifndef TIMER_H
#define TIMER_H

#include <string>

class Timer
{
	public:
		Timer(std::string s) : m_desc(s), m_start(clock()) {}
		~Timer() { 	printf("[Timer] %s %.2lfs\n", m_desc.c_str(), (clock() - m_start) / CLOCKS_PER_SEC); }

	private:
		std::string m_desc;
		double m_start;
};

#endif
