#include "CD_HSE.h"

double CD_BOND::max_static_coord(int dim){
    return std::max(Coords(m_bond->start)[dim],
		    Coords(m_bond->end)[dim]);
}

double CD_BOND::min_static_coord(int dim){
    return std::min(Coords(m_bond->start)[dim],
		    Coords(m_bond->end)[dim]);
}

double CD_BOND::max_moving_coord(int dim,double dt){
    double ans = -HUGE;
    for (int i = 0; i < 2; ++i){
	POINT* pt = (i == 0)? m_bond->start : m_bond->end;
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(ans,sl->x_old[dim]);
	ans = std::max(ans,sl->x_old[dim]+sl->avgVel[dim]*dt); 
    }    
    return ans;
}

double CD_BOND::min_moving_coord(int dim,double dt){
    double ans = HUGE;
    for (int i = 0; i < 2; ++i){
	POINT* pt = (i == 0)? m_bond->start : m_bond->end;
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(ans,sl->x_old[dim]);
	ans = std::min(ans,sl->x_old[dim]+sl->avgVel[dim]*dt); 
    }    
    return ans;
}

POINT* CD_BOND::Point_of_hse(int i) const{
    if (i >= num_pts())
	return NULL;
    else
        return (i == 0) ? m_bond->start : 
			  m_bond->end;
}

double CD_TRI::max_static_coord(int dim){
    double ans = -HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(sl->x_old[dim],ans);
    }
    return ans;
}

double CD_TRI::min_static_coord(int dim){
    double ans = HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(sl->x_old[dim],ans);
    }
    return ans;
}

double CD_TRI::max_moving_coord(int dim,double dt){
    double ans = -HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::max(ans,sl->x_old[dim]);
	ans = std::max(ans,sl->x_old[dim]+sl->avgVel[dim]*dt);
    }
    return ans;
}

double CD_TRI::min_moving_coord(int dim,double dt){
    double ans = HUGE;
    for (int i = 0; i < 3; ++i){
	POINT* pt = Point_of_tri(m_tri)[i];
	STATE* sl = (STATE*)left_state(pt);
	ans = std::min(ans,sl->x_old[dim]);
	ans = std::min(ans,sl->x_old[dim]+sl->avgVel[dim]*dt);
    }
    return ans;
}

POINT* CD_TRI::Point_of_hse(int i) const{
    if (i >= num_pts())
	return NULL;
    else
        return Point_of_tri(m_tri)[i];
}

