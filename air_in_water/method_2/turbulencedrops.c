#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

double tend = 4;
int LEVEL = 7;
double R0 = 0.1;
double We;
double Pe = 1000;

double rrhocp = 1.123 / 4172.4;
double rkappa = 0.026 / 0.6;

scalar T1[]; //temperature 1
scalar T2[]; //temperature 2

double rhocp1 = 1;
double rhocp2;
double kappa1;
double kappa2;

FILE * fpout;

int main (int argc, char *argv[] ) {

  X0 = Y0 = -0.5;
  L0 = 1.;

  size (L0);

  rhocp2 = rhocp1/rrhocp;
  kappa1 = 1/Pe;
  kappa2 = kappa1/rkappa;

  double Re = 1000;
  We = 15;
  double rrho = 1.117 / 997.0;
  double rmu = 0.0000185 / 0.001;

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

  fpout = fopen("output/data.vof", "w");

  int Nvort = 10;

  fraction (f, sq(R0) - (sq(y) + sq(x)) );

  initial_condition(10, u);

  foreach () {
	T1[] = 1.0*f[];
	T2[] = 0.0;
  }

  f.tracers = list_copy ({T1,T2});
  T2.inverse = true;
  
}

event diffusion (i++) {

  scalar T[];
  scalar s[];
  scalar lambda[];
  foreach () {
	double rhocp = rhocp1*f[] + rhocp2*(1. - f[]);
  	T[] = T1[] + T2[];
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
	T1[] = T[]*f[];
	T2[] = T[]*(1.-f[]);
  }

}


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){1.e-2}, LEVEL);
}
#endif


scalar ent[], ke[];
event output (i++) {

  scalar gfrac[], gw[], gw2[];
  foreach () {
    ent[] = rhocp2*T2[];
    ke[] = rho1*f[]*(sq(u.x[]) + sq(u.y[])) + rho2*(1. - f[])*(sq(u.x[]) + sq(u.y[]));
  }


  printf("%g %g %g \n", t, statsf(ent).sum, statsf(ke).sum);
  fprintf(fpout, "%g %g %g \n", t, statsf(ent).sum, statsf(ke).sum);
}



FILE * fp;
event video ( i +=10; t <= tend) {

  vorticity (u, omega);
  output_ppm (f, fp, n=2048, file= "output/f.mp4");
  output_ppm (omega, fp, n=2048, file= "output/omega.mp4");

  scalar T[];
  foreach () {
	double rhocp = rhocp1*f[] + rhocp2*(1. - f[]);
  	T[] = T1[] + T2[];
  }
  output_ppm (T, fp, n=2048, file= "output/T.mp4");

}

