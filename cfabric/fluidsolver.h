#ifndef FLUIDSOLVER_H
#define FLUIDSOLVER_H

struct L_RECTANGLE
{
    int comp {-1};			 
	int m_index {-1};
	double m_coords[MAXD];	
	int icoords[MAXD];

	void setCoords(double* coords, int dim)
    {
        for (int i = 0; i < dim; ++i)
            m_coords[i] = coords[i];
    }
};


class CARTESIAN
{

};


#endif
