#ifndef PROFILER_H
#define PROFILER_H
#include <chrono>
using namespace std;
class StaticProfiler {
    static mutex mtx_;
    static map<int, float> stats_sum;
    static map<int, int> stats_count;
    chrono::steady_clock::time_point startT;
    public:
    void start(){
        startT=chrono::steady_clock::now();
    }
    void elapsed(int i){
        lock_guard<mutex> lk(mtx_);
        float RT = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now()-startT).count();
        stats_sum[i]+=RT;
        stats_count[i]++;
        std::cerr<<"ELAPSED "<<i<<": "<<RT<<"\tM : "<<stats_sum[i]/stats_count[i] <<std::endl;
        startT=chrono::steady_clock::now();
    }
};
class Profiler {
    chrono::steady_clock::time_point startT;
    clock_t begin_time;
    public:
    void start(){
        startT=chrono::steady_clock::now();
        begin_time = clock();
    }
    void check(int i){
        float cpu_millis = float( clock() -begin_time ) / CLOCKS_PER_SEC* 1000;
        float RT = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now()-startT).count();
        float ratio=cpu_millis/float(RT);
        ratio=(8-ratio)*RT;
        std::cerr<<"CHECK "<<i<<": "<<ratio/1000<<std::endl;
        startT=chrono::steady_clock::now();
        begin_time=clock();
    }
    void elapsed(int i){

        float RT = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-startT).count();
        std::cerr<<"ELAPSED "<<i<<": "<<RT<<std::endl;
        startT=chrono::steady_clock::now();
        begin_time=clock();
    }
};
#endif
