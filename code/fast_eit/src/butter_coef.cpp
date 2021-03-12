#include "eit.h"



double butter_coef(unsigned &ii,unsigned &jj, unsigned &kk, double &ximax, double &ximax2,double &steplen,double &steplen2, double &ximaxsteplen, unsigned &N, double omegac)
{
    double v1 = ximax2-2*ii*ximaxsteplen+ii*ii*steplen2;
    double v2 = ximax2-2*jj*ximaxsteplen+jj*jj*steplen2;
    double v3 = ximax2-2*kk*ximaxsteplen+kk*kk*steplen2;

    double omega2 = v1+v2+v3;

    return 1/(sqrt(1+pow(omegac*omega2,N)));
}