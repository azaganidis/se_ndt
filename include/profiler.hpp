#ifndef PROFILER_H
#define PROFILER_H
#include <chrono>
using namespace std;
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

        float RT = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now()-startT).count();
        std::cerr<<"ELAPSED "<<i<<": "<<RT<<std::endl;
        startT=chrono::steady_clock::now();
        begin_time=clock();
    }
};
#endif
