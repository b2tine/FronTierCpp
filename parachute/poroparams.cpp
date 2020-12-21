#include <iostream>
#include <cmath>

double canopyD(double phi, double Dc, double Dp);
double canopyV(double Dp);
double canopyS(double Dc);


int main()
{
    const double Dp = 4.15;       //projected diameter of canopy
    const double Dc = 8.5344;     //structural diameter of canopy

    double h = 0.001;           //canopy thickness
    double rho = 1.29;          //density fluid
    double mu = 0.0001813;      //dynamic viscosity fluid

    
    /*
    double porosity[10] = {0.01,0.02,0.03,0.04,0.05,
                           0.06,0.07,0.08,0.09,0.10};
    */


    double alpha[26];
    double beta[26];


    const double Dsv = Dc*Dc/(Dp*Dp*Dp);
    
    for (int i = 0; i <= 25; ++i)
    {
        double phi = ((double)i)/100.0;
        double phi3 = phi*phi*phi;
        
        alpha[i] = 37.5*mu*Dsv*Dsv/phi3*h;
        beta[i] = 0.875*rho*Dsv/phi3*h;

        /*
        double D = canopyD(phi,Dc,Dp);
        alpha[i] = 150.0*mu*pow((1.0-phi),2.0)/D/D/phi3*h;
        beta[i] = 1.75*rho*(1.0-phi)/D/phi3*h;
        */

        printf("phi = %g:  ",phi);
        printf("alpha = %g , beta = %g\n",alpha[i],beta[i]);
    }

    return 0;
}


//Characteristic length of the canopy
double canopyD(double phi, double Dc, double Dp)
{
    double V = canopyV(Dp);
    double S = canopyS(Dc);
    double D = 6.0*(1.0 - phi)*V/S;
    return D;
}

double canopyV(double Dp)
{
    double V = M_PI*Dp*Dp*Dp/12.0;
    return V;
}

double canopyS(double Dc)
{
    double S = M_PI*Dc*Dc/4.0;
    return S;
}
