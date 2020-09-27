#ifndef CD_HSE_H
#define CD_HSE_H

#include <FronTier.h>
#include "ifluid_state.h"
#include "vtk.h"

#include <string>


enum class CD_HSE_TYPE
{
    FABRIC_TRI,
    STRING_BOND,
    RIGID_TRI
};


//abstract base class for hypersurface element(HSE)
//can be a point, bond, or triangle
struct CD_HSE
{
    CD_HSE_TYPE type;
	virtual double max_static_coord(int) = 0;
	virtual double min_static_coord(int) = 0;
	virtual double max_moving_coord(int,double) = 0;
	virtual double min_moving_coord(int,double) = 0;
	virtual POINT* Point_of_hse(int) const  = 0;
	virtual int num_pts() const= 0;
	virtual ~CD_HSE(){};
};

//wrap class for triangle
struct CD_TRI: public CD_HSE
{
    TRI* m_tri;
	
    CD_TRI(TRI* tri, CD_HSE_TYPE tag)
        : m_tri{tri}
    {
        type = tag;
    }

	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts() const {return 3;}
};

//wrap class for bond
struct CD_BOND: public CD_HSE
{
	int m_dim;
    BOND* m_bond;
	
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

	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const{return 2;}
};

/*
//wrap class for point
struct CD_POINT: public CD_HSE
{
    POINT* m_point;

    CD_POINT(POINT* point)
        : m_point(point)
    {}

	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const {return 1;}
};
*/


bool adjacentHSE(CD_HSE* A, CD_HSE* B);



#endif
