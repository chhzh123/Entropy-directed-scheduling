// Copyright 2018 SYSU
// Author: Hongzheng Chen
// E-mail: chenhzh37@mail2.sysu.edu.cn

// This is the implementation of Entropy-directed scheduling (EDS) algorithm for FPGA high-level synthesis.

// This head file contains the implement of high-accuracy time counter.

#ifndef WATCH_H
#define WATCH_H

// #define NOMINMAX
#include <windows.h>
// #undef NOMINMAX

class stop_watch
{
public:
    stop_watch()
        : elapsed_(0)
    {
        QueryPerformanceFrequency(&freq_);
    }
    ~stop_watch(){}
public:
    void start()
    {
        QueryPerformanceCounter(&begin_time_);
    }
    void stop()
    {
        LARGE_INTEGER end_time;
        QueryPerformanceCounter(&end_time);
        elapsed_ += (end_time.QuadPart - begin_time_.QuadPart) * 1000000 / freq_.QuadPart;
    }
    void restart()
    {
        elapsed_ = 0;
        start();
    }
    //microsecond
    double elapsed()
    {
        return static_cast<double>(elapsed_);
    }
    //millisecond
    double elapsed_ms()
    {
        return elapsed_ / 1000.0;
    }
    //second
    double elapsed_second()
    {
        return elapsed_ / 1000000.0;
    }

private:
    LARGE_INTEGER freq_;
    LARGE_INTEGER begin_time_;
    long long elapsed_;
};

#endif // WATCH_H