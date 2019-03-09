#ifndef FRONT_WRAPPER_H
#define FRONT_WRAPPER_H

#include <FronTier.h>

#include <vector>
#include <string>
#include <iostream>


class FrontWrapper
{
    private:

        Front cfront;

        //Initialization
        void ResetTime();
        void ClipIntfcToSubdomain();


    public:

        //Initialization
        void StartUp(F_BASIC_DATA* ft_basic);
        void InitIntfc(LEVEL_FUNC_PACK* level_func_pack);
        void InitFrontVeloFunc(VELO_FUNC_PACK* velo_func_pack);
        
        //void ReadSpaceDomain(std::string in_name);//,f_basic);
        //void ReadTimeControl(std::string in_name);


        //Manual Setters (move these into a testing only class)
        void setMaxTime(double time);
        void setMaxStep(double step);
        void setPrintTimeInterval(double interval);
        void setMovieFrameInterval(double interval);
        void setFrequencyOfRedistribution(double val, int i = GENERAL_WAVE);


        //Getters
        int Dim();
        double TimeStepFactor();
        //Front* cFront() { return &cfront; }


        //Memory
        void FreeMainIntfc();


        //Mesh/Interface
        void RedistMesh();


        //Output
        bool IsSaveTime();
        bool IsDrawTime();
        void Save();
        void Draw();
        void PrintTimeStamp();

        //Propagation
        void Propagate();
        void SetNextTimeStep();
        void SetOutputCounter();
        void TimeControlFilter();
        void AddTimeStepToCounter();
        bool TimeLimitReached();
};





#endif
