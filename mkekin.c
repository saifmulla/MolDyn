#include <stdio.h>
/*
 *  Scale forces, update velocities and compute K.E.
 */
  double
  mkekin(const int npart, double f[], double vh[], const double hsq2, const double hsq){
    int i;
    double sum=0.0, ekin;

    for (i=0; i<3*npart; i++) {
      f[i]*=hsq2;
      vh[i]+=f[i];
      sum+=vh[i]*vh[i];
    }
    ekin=sum/hsq;

    return(ekin);
  }
