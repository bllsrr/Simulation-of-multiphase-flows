#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

double tend = 5;
int LEVEL = 7;
double R0 = 0.1;
double We;
double Pe = 100;

double rrhocp = 1836.0 / 4172.4;
double rkappa = 0.17 / 0.6;

scalar fh1[]; //enthalpie fluid 1
scalar fh2[]; //enthalpie fluid 2

double rhocp1 = 1;
double rhocp2;
double kappa1;
double kappa2;

FILE * fpout;

int main (int argc, char *argv[] ) {

  if (argc > 1) {
    LEVEL = atoi(argv[1]);
  }

  X0 = Y0 = -0.5;
  L0 = 1.;

  size (L0);

  rhocp2 = rhocp1/rrhocp;
  kappa1 = 1/Pe;
  kappa2 = kappa1/rkappa;

  double Re = 1000;
  We = 5;
  double rrho = 918.0 / 997.0;
  double rmu = 0.084 / 0.001;

  init_grid (1 << LEVEL);
  
  rho1 = 1.0;
  rho2 = rho1/rrho;
  mu1 = 1/Re;
  mu2 = mu1/rmu;
  f.sigma = 1/We;

  /**
  We reduce the tolerance on the Poisson and viscous solvers to
  improve the accuracy. */
  
  TOLERANCE = 1e-4;

  run();

}

scalar omega[];

#include "initialcondition.h"
event init (t = 0) {
  char filename[100];
  sprintf(filename,"output/LEVEL_%d/data.vof",LEVEL);
  fpout = fopen(filename, "w");

  int Nvort = 10;

  fraction (f, sq(R0) - (sq(y) + sq(x)) );

  initial_condition(10, u);

  foreach () {
	fh1[] = rhocp1*f[];
	fh2[] = 0.0;
  }

  f.tracers = list_copy ({fh1,fh2});
  fh2.inverse = true;
  
}

event diffusion (i++) {

  scalar T[];
  scalar s[];
  scalar lambda[];
  foreach () {
	double rhocp = rhocp1*f[] + rhocp2*(1. - f[]);
  	T[] = (fh1[] + fh2[])/rhocp;
  	s[] = -rhocp*T[];
	lambda[] = -rhocp;
  }

  face vector kappa[];
  foreach_face () {
	double ff = (f[1] + f[])/2.;
	kappa.x[] = (kappa1*ff + kappa2*(1. - ff))*dt;
  }

  TOLERANCE = 1.e-4;
  poisson(T, s, kappa, lambda);

  //we reconstruct the enthalpy
  foreach () {
	fh1[] = rhocp1*T[]*f[];
	fh2[] = rhocp2*T[]*(1.-f[]);
  }

}


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){1.e-2}, LEVEL);
}
#endif



scalar T[], T1[], T2[], ke[];
event output (i++) {
  foreach () {
    double rhocp = rhocp1*f[] + rhocp2*(1. - f[]);
    T[] = (fh1[] + fh2[])/rhocp;
  }
  stats sT = statsf(T);
  scalar gfrac[], gw[], gw2[];
  foreach () {
    T1[] = fh1[] / rhocp1*f[];
    T2[] = fh2[] / rhocp2*(1.-f[]);
    ke[] = rho1*f[]*(sq(u.x[]) + sq(u.y[])) + rho2*(1. - f[])*(sq(u.x[]) + sq(u.y[]));
  }

  printf("%g %g %g %g \n", t, sT.max, sT.min, statsf(ke).sum);
  fprintf(fpout, "%g %g %g %g \n", t, sT.max, sT.min, statsf(ke).sum);

  if (sT.max - sT.min < 0.005) {
    return 1; 
  }
}



FILE * fp;
event video (i +=10; t <= tend) {

  vorticity (u, omega);

  char f_mp4[100], omega_mp4[100], T_mp4[100];
  sprintf(f_mp4, "output/LEVEL_%d/f.mp4", LEVEL);
  sprintf(omega_mp4, "output/LEVEL_%d/omega.mp4", LEVEL);
  sprintf(T_mp4, "output/LEVEL_%d/T.mp4", LEVEL);

  output_ppm (f, fp, n=2048, min=0, max=1, file=f_mp4);
  output_ppm (omega, fp, n=2048, min=0, max=1, file=omega_mp4);
  output_ppm (T, fp, n=2048, min=0, max=1, file=T_mp4);

}