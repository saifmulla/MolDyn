#include "MolDyn.h"
/*
 *  Variable declarations
 */

    double epot;
    double vir;
    double count;

/*
 *  Main program : Molecular Dynamics simulation.
 */
int main(int argc, char* argv[]){
    //initial sizes
    int mm = 12;// size of box
    int i = 0;
    double h = 0.0025;
    int movemx = 20;
    
    /**
     * check if command line values are supplied if yes
     * then use those values else
     * stay with default values
     */
    if (argc > 1) {
        if(argv[1])
            mm = atoi(argv[1]);
        //check time step size
        if (argv[2])
            h = (double) atof(argv[2]);
        if (argv[3])
            movemx = atoi(argv[3]);
    }
    
    int npart = 4*mm*mm*mm; //size of molecules
    int move;
    double* x=NULL;
    double* vh=NULL;
    double* f=NULL;
    x = (double*) malloc(sizeof(double)*npart*3);
    vh = (double*) malloc(sizeof(double)*npart*3);
    f = (double*) malloc(sizeof(double)*npart*3);
    
    unsigned long memalign16 = ((unsigned long)x & 15);
    unsigned long memalign32 = ((unsigned long)vh & 31);
    
    printf("Memalign16 %s\n",(memalign16==0?"aligned":"unaligned"));
    printf("Memalign32 %s\n",(memalign32==0?"aligned":"unaligned"));
    
    double ekin;
    double vel;
    double sc;
    double start, time;


  /*
   *  Parameter definitions
   */

    double den    = 0.83134;
    double side   = pow((double)npart/den,0.3333333);
    double tref   = 0.722;
    double rcoff  = (double)mm/4.0;
    int    irep   = 10;
    int    istop  = 20;
    int    iprint = 5;

    double a      = side/(double)mm;
    double hsq    = h*h;
    double hsq2   = hsq*0.5;
    double tscale = 16.0/((double)npart-1.0);
    double vaver  = 1.13*sqrt(tref/24.0);

  /*
   *  Initial output
   */

    printf(" Molecular Dynamics Simulation example program\n");
    printf(" ---------------------------------------------\n");
    printf(" number of particles is ............ %6d\n",npart);
    printf(" side length of the box is ......... %13.6f\n",side);
    printf(" cut off is ........................ %13.6f\n",rcoff);
    printf(" reduced temperature is ............ %13.6f\n",tref);
    printf(" basic timestep is ................. %13.6f\n",h);
    printf(" temperature scale interval ........ %6d\n",irep);
    printf(" stop scaling at move .............. %6d\n",istop);
    printf(" print interval .................... %6d\n",iprint);
    printf(" total no. of steps ................ %6d\n",movemx);

  /*
   *  Generate fcc lattice for atoms inside box
   */
    fcc(x, npart, mm, a);
  /*
   *  Initialise velocities and forces (which are zero in fcc positions)
   */
    mxwell(vh, 3*npart, h, tref);
    dfill(3*npart, 0.0, f, 1);
  /*
   *  Start of md
   */
    printf("\n    i       ke         pe            e         temp   "
           "   pres      vel      rp\n  -----  ----------  ----------"
           "  ----------  --------  --------  --------  ----\n");

     start = secnds(); 


    for (move=1; move<=movemx; move++) {

    /*
     *  Move the particles and partially update velocities
     */
      domove(3*npart, x, vh, f, side);

    /*
     *  Compute forces in the new positions and accumulate the virial
     *  and potential energy.
     */
      forces(npart, x, f, side, rcoff);

    /*
     *  Scale forces, complete update of velocities and compute k.e.
     */
      ekin=mkekin(npart, f, vh, hsq2, hsq);
    /*
     *  Average the velocity and temperature scale if desired
     */
      vel=velavg(npart, vh, vaver, h);
      if (move<istop && fmod(move, irep)==0) {
        sc=sqrt(tref/(tscale*ekin));
        dscal(3*npart, sc, vh, 1);
        ekin=tref/tscale;
      }

    /*
     *  Sum to get full potential energy and virial
     */
      if (fmod(move, iprint)==0)
        prnout(move, ekin, epot, tscale, vir, vel, count, npart, den);
      
    }

    time = secnds() - start;  
#ifdef __openmp
    printf("Time =  %f\n",(float) time);
#else
    printf("Time = %f\n",((float) time)/CLOCKS_PER_SEC);
#endif
  }

time_t starttime = 0; 

double secnds()
{
#ifdef __openmp
  return omp_get_wtime();
#else
    return clock();
#endif

}

