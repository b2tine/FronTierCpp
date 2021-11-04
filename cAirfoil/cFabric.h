#ifndef CFABRIC_H
#define CFABRIC_H

#include <cFluid.h>


class CFABRIC_CARTESIAN : public G_CARTESIAN
{
public:

    CFABRIC_CARTESIAN(Front* ft)
        : G_CARTESIAN(ft)
    {}

    ~CFABRIC_CARTESIAN() = default;

    void applicationSetComponent();

    void applicationSetStates();

protected:

    void addFluxAlongGridLine(int,int*,double,SWEEP*,FSWEEP*) override;

    void appendGhostBuffer(SWEEP*,SWEEP*,int,int*,int,int) override;

private:

    void setElasticStates(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                        int*,int,int,int,int,int);

    void setElasticStatesRFB_normal(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

    void setElasticStatesRFB(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

    void setElasticStatesRiem(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);
};


#endif
