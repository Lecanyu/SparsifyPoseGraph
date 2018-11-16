/* MIT License
 *
 * Copyright (c) 2018 Björn Blissing
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef USE_PROFILER
#define USE_PROFILER
#endif

#ifndef _FLAMEPROFILER_H_
#define _FLAMEPROFILER_H_

#ifndef USE_PROFILER
// Expand to nothing if profiling is disabled
#define PZone(name)
#define PZoneCat( name, category )
#define PMetadata( title, value )
#define ProfilerWriteFileHeadToDisk()
#else
// Macro concatenation
#define CONCAT_IMPL(x, y) x##y
#define MACRO_CONCAT(x, y) CONCAT_IMPL( x, y )
// Macro names
#define PZone(name) Profiler::Zone MACRO_CONCAT( profileZone, __COUNTER__ )( name )
/**
 * Canyu Le
 * Allow user to set profile output path
 * */
// -----------------------
#define PZonePath(name, category, path, filename) Profiler::Zone MACRO_CONCAT( profileZone, __COUNTER__ )( name, category, path, filename )
// -----------------------
#define PZoneCat(name, category) Profiler::Zone MACRO_CONCAT( profileZone, __COUNTER__ )( name , category )
#define PMetadata(title, value) Profiler::FlameGraphWriter::instance().addMetadata( title, value )
#define ProfilerWriteFileHeadToDisk() Profiler::FlameGraphWriter::instance().WriteBufferedDataToDisk()
#define ProfilerWriteFileHeadToDiskWithPath(path, filename) Profiler::FlameGraphWriter::instance(path, filename).WriteBufferedDataToDisk()

#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <algorithm>
#include <limits>
#include <thread>
#include <cstdio>
#include <chrono>
#include <iostream>
#include <atomic>

#ifdef _WIN32
#include <windows.h>
#else

#include <unistd.h>

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sys/time.h>
#include <sys/resource.h>
#include <boost/chrono/duration.hpp>
#include <boost/chrono/thread_clock.hpp>
#include <boost/chrono.hpp>

#endif

namespace Profiler
{

    struct TracePoint
    {
        char name[64];
        char category[40];
        uint64_t timeStart;
        uint64_t timeEnd;
        uint32_t processId;
        uint32_t threadId;
        int32_t ru_nivcsw;
    }; // 128 byte struct

    class FlameGraphWriter
    {
    public:

        ~FlameGraphWriter()
        {
            std::cout << "~FlameGraphWriter()" << std::endl;
            WriteToDiskBulk();
        }

        // for some reason, the de-constructor is not called on robot, so we define a explict
        void WriteToDiskBulk()
        {
            if(file_tail_writed)
            {
                file_writed = true;
            }

            if (file_writed)
            {
                return;
            }
            file_writed = true;

            if(!file_head_writed)
            {
                WriteFileHeadToDisk();
                file_head_writed = true;
            }
            WriteBufferedDataToDisk();
            if(!file_tail_writed)
            {
                WriteFileTailToDisk();
                file_tail_writed = true;
            }

        }

        void WriteFileHeadToDisk()
        {
            if(file_head_writed)
            {
                return;
            }

            if(!ofs.is_open())
            {
                std::string file_full_name = file_path + "/" + file_name;
                // Open file
                ofs.open(file_full_name, std::ofstream::out);
                std::cout << "====================================" << std::endl;
            }


            // Header
            ofs << "{\n";

            for (const auto &metadata : m_metadata)
            {
                ofs << "\t\"" << metadata.first << "\": \"" << metadata.second << "\"";
                ofs << ",\n";
            }

            ofs << "\"traceEvents\": ";
            // Begin trace events
            ofs << "[";

            ofs.flush();
            file_head_writed = true;
        }

        void WriteFileTailToDisk()
        {
            if(file_tail_writed)
            {
                return;
            }

            const auto &tracepoint = last_trace_point;
            {
                uint64_t duration = tracepoint.timeEnd - tracepoint.timeStart;
                ofs << "\n";
                ofs << "{";
                ofs << " \"pid\":" << tracepoint.processId << ",";
                ofs << " \"tid\":" << tracepoint.threadId << ",";
                ofs << " \"ts\":" << tracepoint.timeStart - profiler_start_time << ",";
                ofs << " \"dur\":" << duration << ",";
                ofs << " \"ph\":\"X\",";
                ofs << " \"name\":\"" << tracepoint.name << "\",";
                ofs << " \"cat\":\"" << tracepoint.category << "\",";

                ofs << " \"args\": ";
                ofs << "\t{";
                ofs << "\"ru_nivcsw\":";
                ofs << tracepoint.ru_nivcsw ;

                ofs << "}";

                ofs << "}";
            }

            ofs << "\n\t]";
            ofs << "\n}\n";
            file_tail_writed = true;
            ofs.flush();
            ofs.close();
        }

        void WriteBufferedDataToDisk()
        {
            if (!file_head_writed)
            {
                WriteFileHeadToDisk();
                file_head_writed = true;
            }

            for (const auto &tracepoint : m_tracepoints)
            {
                profiler_start_time = (std::min)(profiler_start_time, tracepoint.timeStart);
            }

            // Write content
            if (ofs.is_open())
            {
                // For each tracepoint write

                for (const auto &tracepoint : m_tracepoints)
                {
                    uint64_t duration = tracepoint.timeEnd - tracepoint.timeStart;
                    ofs << "\n";
                    ofs << "{";
                    ofs << " \"pid\":" << tracepoint.processId << ",";
                    ofs << " \"tid\":" << tracepoint.threadId << ",";
                    ofs << " \"ts\":" << tracepoint.timeStart - profiler_start_time << ",";
                    ofs << " \"dur\":" << duration << ",";
                    ofs << " \"ph\":\"X\",";
                    ofs << " \"name\":\"" << tracepoint.name << "\",";
                    ofs << " \"cat\":\"" << tracepoint.category << "\",";

                    ofs << " \"args\": ";
                    ofs << "\t{";
                    ofs << "\"ru_nivcsw\":";
                    ofs << tracepoint.ru_nivcsw ;

                    ofs << "}";


                    ofs << "}";
                    ofs << ",";
                }

                if(!m_tracepoints.empty())
                {
                    last_trace_point =m_tracepoints.back();
                    m_tracepoints.clear();
                }

                ofs.flush();
            }
        }

        static FlameGraphWriter &instance()
        {
            static FlameGraphWriter instance;
            return instance;
        }

        static FlameGraphWriter &instance(std::string dir, std::string filename)
        {
            static FlameGraphWriter instance(dir, filename);
            return instance;
        }

        void addMetadata(const std::string &title, const std::string &value)
        {
            std::pair<std::string, std::string> metadata = std::make_pair(title, value);
            m_metadata.push_back(metadata);
        }

        void addTracePoint(const TracePoint &point)
        {
            m_tracepoints.push_back(point);
        }

        FlameGraphWriter(const FlameGraphWriter &) = delete; // No copy allowed
        FlameGraphWriter(const FlameGraphWriter &&) = delete; // No move allowed
        FlameGraphWriter &operator=(const FlameGraphWriter &) = delete; // No assignment allowed
        FlameGraphWriter &operator=(FlameGraphWriter &&) = delete; // No move assignment allowed

        FlameGraphWriter(const std::string& dir, const std::string& filename)
        {
            file_path = dir;
            file_name = filename;

            // Reserve points to avoid early memory reallocations
            m_tracepoints.reserve(1000);

            file_writed = false;
            file_head_writed = false;
            file_tail_writed = false;
        }
    private:
        std::ofstream ofs;
        std::string file_path = "/home/lecanyu/Desktop/";
        std::string file_name = "graph_localization_profiler.json";

        FlameGraphWriter()
        {
            // Reserve points to avoid early memory reallocations
            m_tracepoints.reserve(1000);

            file_writed = false;
            file_head_writed = false;
            file_tail_writed = false;
        }

        std::vector<TracePoint> m_tracepoints;
        std::vector<std::pair<std::string, std::string>> m_metadata;

        std::atomic_bool file_writed;
        std::atomic_bool file_head_writed;
        std::atomic_bool file_tail_writed;
        uint64_t profiler_start_time = (std::numeric_limits<uint64_t>::max)();
        TracePoint last_trace_point;
    };

    class Zone
    {
    public:
        /**
         * Canyu Le
         * Allow user to set profile output path
         * */
        // -----------------------
        Zone(const char *name, const char *category = "default", const char* dir = "None", const char* filename = "None")
            :_dir(dir), _filename(filename)
        {
            getrusage(RUSAGE_SELF, &usage_start);

            start_time=boost::chrono::thread_clock::now();

            snprintf(m_tracepoint.name, sizeof(m_tracepoint.name), "%s", name);
            snprintf(m_tracepoint.category, sizeof(m_tracepoint.category), "%s", category);
            m_tracepoint.threadId = 1;
            auto timer = std::chrono::high_resolution_clock::now();
            m_tracepoint.timeStart = std::chrono::duration_cast<std::chrono::microseconds>(
                    timer.time_since_epoch()).count();
            m_tracepoint.processId = static_cast<uint32_t>(::getpid());

            // show on cpu time in same time-line start with same start time
            snprintf(m_tracepoint_on_cpu.name, sizeof(m_tracepoint_on_cpu.name), "%s", name);
            snprintf(m_tracepoint_on_cpu.category, sizeof(m_tracepoint_on_cpu.category), "%s", category);
            m_tracepoint_on_cpu.threadId = m_tracepoint.threadId + 1;
            m_tracepoint_on_cpu.timeStart = m_tracepoint.timeStart;
            m_tracepoint_on_cpu.processId = m_tracepoint.processId + 1;
        }
        // -----------------------

        ~Zone()
        {
            auto timer = std::chrono::high_resolution_clock::now();
            m_tracepoint.timeEnd = std::chrono::duration_cast<std::chrono::microseconds>(
                    timer.time_since_epoch()).count();
            {
                struct rusage usage_end;
                getrusage(RUSAGE_SELF, &usage_end); // precision: 1 us on PC, 5~10 ms on board.

                boost::chrono::thread_clock::duration elapsed = boost::chrono::thread_clock::now() - start_time;
                uint64_t dt = boost::chrono::nanoseconds(elapsed).count()/1e3;

                m_tracepoint_on_cpu.timeEnd = m_tracepoint.timeStart + dt;

                m_tracepoint_on_cpu.ru_nivcsw = usage_end.ru_nivcsw - usage_start.ru_nivcsw;

                m_tracepoint.ru_nivcsw = m_tracepoint_on_cpu.ru_nivcsw;
                /**
                 * Canyu Le
                 * Allow user to set profile output path
                 * */
                // -----------------------
                if(_dir=="None" && _filename=="None")
                {
                    FlameGraphWriter::instance().addTracePoint(m_tracepoint);
                    FlameGraphWriter::instance().addTracePoint(m_tracepoint_on_cpu);
                }
                else
                {
                    FlameGraphWriter::instance(_dir, _filename).addTracePoint(m_tracepoint);
                    FlameGraphWriter::instance(_dir, _filename).addTracePoint(m_tracepoint_on_cpu);
                }
                // -----------------------
            }
        }

    private:
        TracePoint m_tracepoint;

        // add trace point that just account on CPU time
        TracePoint m_tracepoint_on_cpu;
        struct rusage usage_start;

        boost::chrono::thread_clock::time_point start_time; // counting thread on-CPU time.

        // path
        std::string _dir, _filename;
    };
} // namespace profiler

#endif // USE_PROFILER

#endif // !_FLAMEPROFILER_H_
