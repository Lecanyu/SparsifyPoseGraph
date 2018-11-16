/*
 * stopwatch.h
 *
 *  Created on: 23/mar/2013
 *      Author: Mladen Mazuran
 */

#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#ifdef __APPLE__
#   include <mach/mach_time.h>
#elif defined(__linux)
#   include <time.h>
#elif defined(_WIN32) || defined(_WIN64)
#   include <windows.h>
#endif

#ifdef __APPLE__
class Stopwatch {
    uint64_t total, tstart;
    double conversion;

public:
    inline Stopwatch() : total(0) {
        mach_timebase_info_data_t info;
        mach_timebase_info(&info);
        conversion = info.numer * 1e-9 / info.denom;
    }

    inline Stopwatch(const Stopwatch &s) :
        total(s.total), tstart(s.tstart), conversion(s.conversion) {
    }

    inline Stopwatch &operator=(const Stopwatch &s) {
        total = s.total;
        tstart = s.tstart;
        conversion = s.conversion;
        return *this;
    }

    inline void start() {
        tstart = mach_absolute_time();
    }

    inline void stop() {
        uint64_t tend = mach_absolute_time();
        total += tend - tstart;
    }

    inline void reset() {
        total = 0;
    }

    inline double time() {
        return conversion * total;
    }

};
#elif defined(__linux)
class Stopwatch {
    struct timespec tstart, tend;
    double total;

public:
    inline Stopwatch() : total(0) {
    }

    inline Stopwatch(const Stopwatch &s) :
        tstart(s.tstart), total(total) {
    }

    inline Stopwatch &operator=(const Stopwatch &s) {
        tstart = s.tstart;
        total = s.total;
        return *this;
    }

    inline void start() {
        clock_gettime(CLOCK_MONOTONIC, &tstart);
    }

    inline void stop() {
        clock_gettime(CLOCK_MONOTONIC, &tend);
        total += diff();
    }

    inline void reset() {
        total = 0;
    }

    inline double time() {
        return total;
    }

private:
    double diff()
    {
        double te = (double) tend.tv_sec + (double) tend.tv_nsec * 1e-9;
        double ts = (double) tstart.tv_sec + (double) tstart.tv_nsec * 1e-9;
        return te - ts;
    }

};
#elif defined(_WIN32) || defined(_WIN64)
class Stopwatch {
    __int64 total, tstart;
    double conversion;

public:
    inline Stopwatch() : total(0) {
        LARGE_INTEGER li;
        QueryPerformanceFrequency(&li);
        conversion = 1. / li.QuadPart;
    }

    inline Stopwatch(const Stopwatch &s) :
        total(s.total), tstart(s.tstart), conversion(s.conversion) {
    }

    inline Stopwatch &operator=(const Stopwatch &s) {
        total = s.total;
        tstart = s.tstart;
        conversion = s.conversion;
        return *this;
    }
    inline void start() {
        LARGE_INTEGER li;
        QueryPerformanceCounter(&li);
        tstart = li.QuadPart;
    }

    inline void stop() {
        LARGE_INTEGER li;
        QueryPerformanceCounter(&li);
        total += li.QuadPart - tstart;
    }

    inline void reset() {
        total = 0;
    }

    inline double time() {
        return conversion * total;
    }

};
#else
#   error Missing implementation of Stopwatch for this platform
#endif

#endif /* STOPWATCH_H_ */
