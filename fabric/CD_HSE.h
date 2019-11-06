#ifndef CD_HSE_H
#define CD_HSE_H

#include <FronTier.h>
#include "state.h"

#include <vector>
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
	virtual const POINT* const Point_of_hse(int) const  = 0;
	virtual int num_pts() const= 0;
	virtual ~CD_HSE() = default;
};

//wrap class for triangle
struct CD_TRI: public CD_HSE
{
    TRI* m_tri;
	
    CD_TRI(TRI* tri, const char* n)
        : m_tri{tri}
    {
        name = n;
    }

	double max_static_coord(int) override;
	double min_static_coord(int) override;
	double max_moving_coord(int,double) override;
	double min_moving_coord(int,double) override;
	const POINT* const Point_of_hse(int) const override;
	int num_pts() const noexcept override {return 3;}
};

//wrap class for bond
struct CD_BOND: public CD_HSE
{
    BOND* m_bond;
	
    CD_BOND(BOND* bond, const char* n)
        : m_bond{bond}
    {
        name = n;
    }

	double max_static_coord(int) override;
	double min_static_coord(int) override;
	double max_moving_coord(int,double) override;
	double min_moving_coord(int,double) override;
	const POINT* const Point_of_hse(int) const override;
	int num_pts() const noexcept override {return 2;}
};


const bool AreAdjacentHSEs(const CD_HSE* const A, const CD_HSE* const B);


#endif
