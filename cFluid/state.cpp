#include "state.h"


extern double getStateDens(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->dens;
}	/* end getStateDens */

extern double getStateEngy(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->engy;
}	/* end getStateEngy */

extern double getStatePres(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->pres;
}	/* end getStatePres */

extern double getStateMu(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->mu;
}	/* end getStateMu */

extern double getStateKTurb(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->k_turb;
}	/* end getStateKTurb */

extern double getStateTemp(POINTER state)
{
    STATE *fstate = (STATE*)state;
    return fstate->temp;
}   /* end getStateTemp */

extern double getStateXmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[0];
}	/* end getStateXmom */

extern double getStateYmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[1];
}	/* end getStateYmom */

extern double getStateZmom(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->momn[2];
}	/* end getStateZmom */

extern double getStateXvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[0];
}	/* end getStateXvel */

extern double getStateYvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[1];
}	/* end getStateYvel */

extern double getStateZvel(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vel[2];
}	/* end getStateZvel */

//TODO: Need 3d version also
extern double getStateVort(POINTER state)
{
	STATE *fstate = (STATE*)state;
	return fstate->vort;
}	/* end getStateVort */

extern double getStateXimp(POINTER state)
{
    STATE *fstate = (STATE*)state;
    return fstate->impulse[0];
}   /* end getStateXimp */

extern double getStateYimp(POINTER state)
{
    STATE *fstate = (STATE*)state;
    return fstate->impulse[1];
}   /* end getStateYimp */

extern double getStateZimp(POINTER state)
{
    STATE *fstate = (STATE*)state;
    return fstate->impulse[2];
}   /* end getStateZimp */

