//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Timer - Simple object to that prints the length of its lifetime
//
#include "Timer.h"

Timer::~Timer()
{
	printf("%s -- %.2lf s\n", m_desc.c_str(), (clock() - m_start) / CLOCKS_PER_SEC);
}
