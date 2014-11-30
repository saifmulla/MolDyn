//
//  MolDyn.h
//  MolDyn
//
//  Created by Saif Mulla on 27/11/2014.
//  Copyright (c) 2014 Saif Mulla. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
//check if compiler is GCC
// coz for the moment clang does not contain openmp
#ifdef __openmp
#include <omp.h>
#endif

/*
 *  Function declarations
 */
#ifndef MolDyn_h
#define MolDyn_h


void
dfill(const int,const double,double[],const int);

void
domove(const int,double[],double[],double[],const double);

void
dscal(int,double,double[],int);

void
fcc(double[],const int,const int,const double);

void
forces(const int,double[],double[],const double, const double);

double
mkekin(const int,double[],double[],const double,const double);

void
mxwell(double[],int,double,double);

void
prnout(int,double,double,double,double,double,double,int,double);

double
velavg(const int,double[],const double,const double);

double
secnds(void);

#endif