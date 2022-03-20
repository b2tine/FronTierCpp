#ifndef CFABRIC_H
#define CFABRIC_H

#include <cFluid.h>


class CFABRIC_CARTESIAN : public G_CARTESIAN
{
public:

    explicit CFABRIC_CARTESIAN(Front* ft)
        : G_CARTESIAN{ft}
    {}

    ~CFABRIC_CARTESIAN() = default;

    void applicationSetComponent();

    void applicationSetStates();

protected:

    void addFluxAlongGridLine(int,int*,double,SWEEP*,FSWEEP*) override;

    void appendGhostBuffer(SWEEP*,SWEEP*,int,int*,int,int) override;

    void setViscousGhostState(int* icoords, COMPONENT comp, VSWEEP* vs, SWEEP* m_vst) override;

private:

    void setElasticStates(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                        int*,int,int,int,int,int);

    void setElasticStatesDarcy(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

    void setElasticStatesRFB_normal(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

    void setElasticStatesRFB(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

    void setElasticStatesRiem(SWEEP*,SWEEP*,HYPER_SURF*,STATE*,
                            int*,int,int,int,int,int);

    void setElasticViscousGhostState(int* icoords, SWEEP* m_vst, VSWEEP* vs, double* ghost_coords,
            double* crx_coords, COMPONENT comp, double* intrp_coeffs,
            HYPER_SURF_ELEMENT* hse, HYPER_SURF* hs);
};


#endif
