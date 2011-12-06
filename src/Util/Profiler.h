///----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// Profiler.h -- Lightweight macro-based function profiler.
// Currently not thread safe
//
#ifndef PROFILER_H
#define PROFILER_H

#include <time.h>
#include <iostream>

//#define USE_PROFILER 1
#if defined(HAVE_CLOCK_GETTIME) && defined(USE_PROFILER)

// This class writes the lifespan of the object
// to the output variable, in nanoseconds
class TimeTracker
{
    public:
        TimeTracker(size_t & output) : m_output(output)
        {
             timespec start;
             clock_gettime(CLOCK_REALTIME, &start);
             m_start_ns = start.tv_sec * 1000000000 + start.tv_nsec;
        }

        ~TimeTracker()
        {
             timespec end;
             clock_gettime(CLOCK_REALTIME, &end);
             size_t end_ns = end.tv_sec * 1000000000 + end.tv_nsec;

             // Update the result
             m_output += end_ns - m_start_ns;
        }

    private:
        size_t m_start_ns;
        size_t& m_output;
};

// Place the following macros before/after the code in which you want to profile
// Example:
// PROFILE_START("my_function")
// ... some code
// PROFILE_END
#define PROFILE_TICKS_BEFORE_PRINT 100

#define PROFILE_FUNC(x) static std::string __profile_name = x; \
                         static size_t __profile_iterations = 0; \
                         static size_t __profile_total_nanoseconds = 0; \
                         double milli_seconds = (double)__profile_total_nanoseconds / 1000000.0f; \
                         double avg_per_iteration = milli_seconds / __profile_iterations; \
                         if(__profile_iterations++ % PROFILE_TICKS_BEFORE_PRINT == 0) \
                             printf("[Profile %s] count: %zu time: %.0lf ms avg: %.0lf ms\n", __profile_name.c_str(), __profile_iterations, milli_seconds, avg_per_iteration); \
                         TimeTracker __profile_timer(__profile_total_nanoseconds);
                         
#else

#define PROFILE_FUNC(x)

#endif // #ifdef HAVE_CLOCK_GETTIME

#endif // #ifndef PROFILER_H
