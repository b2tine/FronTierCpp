#ifndef CD_HSE_H
#define CD_HSE_H

#include <FronTier.h>
#include "state.h"

#include <string>


//abstract base class for hypersurface element(HSE)
//can be a point, bond, or triangle
struct CD_HSE
{
    std::string name;
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
	
    CD_TRI(TRI* tri, const char* n)
        : m_tri(tri)
    {
        name = n;
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
	
    CD_BOND(BOND* bond, int dim, const char* n)
        : m_bond(bond), m_dim(dim)
    {
        name = n;
    }

	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const{return 2;}
};

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



bool adjacentHSE(CD_HSE* A, CD_HSE* B);



#endif
