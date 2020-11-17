#ifndef INJECTION_H
#define INJECTION_H

#include <limits>
#include <queue>

class IF_Injection
{
    protected:
    
        double t0;

    public:

        IF_Injection(double t) : t0{t} {}
        const double start() const {return t0;}

        IF_Injection() = default;
        IF_Injection(const IF_Injection&) = default;
        IF_Injection& operator=(const IF_Injection&) = default; 
        IF_Injection(IF_Injection&&) = default;
        IF_Injection& operator=(IF_Injection&&) = default; 
        virtual ~IF_Injection() = default;
};


/*
struct IF_Injection_Less
{
    bool operator(T* a, T* b)
    {
        return (a->start() < b->start()) ? a : b;
    }
};

*/


class VortexInjection
    : public IF_Injection
{
    private:

        double center[MAXD];            // center of vortex
        double D;                       // size of vortex
        double A;                       // intensity of vortex

    public:

//        void operator(double* coords)
};



#endif
