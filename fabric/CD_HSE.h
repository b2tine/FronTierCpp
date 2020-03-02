#ifndef CD_HSE_H
#define CD_HSE_H

#include <FronTier.h>
#include "state.h"

#include <string>


enum class CD_HSE_TYPE
{
    FABRIC_TRI,
    STRING_BOND,
    RIGID_TRI
};

struct CD_HSE
{
    CD_HSE_TYPE type;
	virtual double max_static_coord(int dim) const = 0;
	virtual double min_static_coord(int dim) const = 0;
	virtual double max_moving_coord(int dim, double dt) const = 0;
	virtual double min_moving_coord(int dim, double dt) const = 0;
	virtual POINT* Point_of_hse(int i) const  = 0;
	virtual int num_pts() const noexcept = 0;
	virtual ~CD_HSE(){};
};

struct CD_BOND: public CD_HSE
{
	int m_dim; //remove
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

	double max_static_coord(int dim) const override;
	double min_static_coord(int dim) const override;
	double max_moving_coord(int dim, double dt) const override;
	double min_moving_coord(int dim, double dt) const override;
	POINT* Point_of_hse(int i) const override;
	int num_pts() const noexcept override {return 2;}
};

struct CD_TRI: public CD_HSE
{
    TRI* m_tri;
	
    CD_TRI(TRI* tri, CD_HSE_TYPE tag)
        : m_tri{tri}
    {
        type = tag;
    }

	double max_static_coord(int dim) const override;
	double min_static_coord(int dim) const override;
	double max_moving_coord(int dim, double dt) const override;
	double min_moving_coord(int dim, double dt) const override;
	POINT* Point_of_hse(int i) const override;
	int num_pts() const noexcept override {return 3;}
};


bool adjacentHSE(CD_HSE* A, CD_HSE* B);


#endif
