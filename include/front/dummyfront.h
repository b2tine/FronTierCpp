#include <fdecs.h>

#ifndef DUMMY_FRONT_H
#define DUMMY_FRONT_H

#include <vector>
#include <string>

namespace cpp
{

constexpr int MAXD {3};

struct BasicData
{
    bool ReadFromInput;
    std::string in_name;
    std::string out_name;

    bool ResetTime;
    bool RestartRun;
    int RestartStep;
    std::string restart_name;
    std::string restart_state_name;

    double L[MAXD];
    double U[MAXD];
    int gmax[MAXD];
    int subdomains[MAXD];
    int boundary[MAXD][2];
    size_t size_of_intfc_state;

    GEOMETRY_REMAP coord_system;
};


class Front
{
    public:

        void testfunction();

};


}



#endif
