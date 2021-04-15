#ifndef CD_HSE_H
#define CD_HSE_H

#include <FronTier.h>
#include "state.h"
#include "vtk.h"

#include <string>
#include <utility>
#include <numeric>

//TODO: ideally we should combine SPRING_VERTEX and CD_POINT into
//      a single class that the spring solver and collision solver
//      can both operate on


//Proxy class for POINT
struct CD_POINT
{
    POINT* pt;

    /*
    explicit CD_POINT(POINT* point)
        : pt(point)
    {}
    */

    double collsn_dt;
    double collsnImpulse[3];
    double collsnImpulse_RG[3];
    double strainImpulse[3];
    double friction[3];
    double avgVel[3];
    double avgVel_old[3];
    double x_old[3];
    int strain_num;
    int collsn_num;
    int collsn_num_RG;
    bool has_collsn;
    bool has_strainlim;
    bool is_fixed;
    bool is_movableRG;
    bool is_stringpt;

	struct UF   
    {
        int num_pts;
        CD_POINT* root;
        CD_POINT* tail;
        CD_POINT* next_pt;
    };
    
    UF impZone;
};


enum class CD_HSE_TYPE
{
    FABRIC_TRI,
    STRING_BOND,
    STATIC_RIGID_TRI,
    MOVABLE_RIGID_TRI
};


//abstract base class for hypersurface element(HSE)
//can be a point, bond, or triangle
struct CD_HSE
{
    CD_HSE_TYPE type;
	
    virtual ~CD_HSE() {};

    virtual int num_pts() const= 0;
    virtual double max_static_coord(int) = 0;
	virtual double min_static_coord(int) = 0;
	virtual double max_moving_coord(int,double) = 0;
	virtual double min_moving_coord(int,double) = 0;
	virtual POINT* Point_of_hse(int) const = 0;
	virtual CD_POINT* CD_Point_of_hse(int) = 0;
};

//wrap class for triangle
struct CD_TRI: public CD_HSE
{
    TRI* m_tri;
    CD_POINT cpts[3];

    CD_TRI(TRI* tri, CD_HSE_TYPE tag)
        : m_tri{tri}
    {
        type = tag;
    }

	int num_pts() const {return 3;}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	CD_POINT* CD_Point_of_hse(int);
};

//wrap class for bond
struct CD_BOND: public CD_HSE
{
    BOND* m_bond;
    CD_POINT cpts[2];
	
    CD_BOND(BOND* bond, CD_HSE_TYPE tag)
        : m_bond{bond}
    {
        type = tag;
        if (type == CD_HSE_TYPE::STRING_BOND)
        {
            STATE* sl;
            sl = (STATE*)left_state(m_bond->start);
            sl->is_stringpt = true;
            sl = (STATE*)left_state(m_bond->end);
            sl->is_stringpt = true;
        }
    }

	int num_pts() const{return 2;}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	CD_POINT* CD_Point_of_hse(int);
};


bool adjacentHSE(CD_HSE* A, CD_HSE* B);



#endif
