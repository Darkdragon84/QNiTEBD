#ifndef TICTOC_H
#define TICTOC_H

#include <sys/time.h>
#include <ctime>

//using namespace std;
class tictoc
{
    public:
        tictoc():running_(false),etime_(0) {}
        void tic()
        {
            if(running_) std::cerr<<"already running"<<std::endl;
            else
            {
                etime_=0;
                gettimeofday(&tvs_,NULL);
                running_=true;
            }
        }
        double toc()
        {
            if (running_)
            {
                gettimeofday(&tve_,NULL);
                etime_=(tve_.tv_sec - tvs_.tv_sec) + (double)(tve_.tv_usec - tvs_.tv_usec)/1e6;
                running_=false;
                return etime_;
            }
            else
            {
                std::cerr<<"not running!"<<std::endl;
                return 0;
            }
        }
        virtual ~tictoc() {}
    protected:
    private:
        bool running_;
        double etime_;
        struct timeval tvs_,tve_;
};

#endif // TICTOC_H
