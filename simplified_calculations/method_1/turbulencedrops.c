#include "grid/multigrid.h"
#include "run.h"
#include "fractions.h"
#include "poisson.h"
#include "redistance.h"
#include <stdlib.h>


double tend = 5;
int LEVEL;
double R0 = 0.1;

double rrhocp = 1836.0 / 4172.4;
double rkappa = 0.17 / 0.6;

double Delta1;
int current_m;

scalar d[];
scalar sf[];
scalar fh1[];
scalar fh2[];

double rhocp1 = 1;
double rhocp2;
double kappa1;
double kappa2;

FILE * fpout = NULL;

int main (int argc, char * argv[]) {
  
  int m_val = 4;
  int n_val = 8;

  if (argc == 3){
    m_val = atoi(argv[1]);
    n_val = atoi(argv[2]);
  } 
  
  else {
    fprintf(stderr, "Utilisation: ./a.out <m> <n>\n");
    return 1;
  }

  X0 = Y0 = -0.5;
  L0 = 1.;

  size (L0);

  rhocp2 = rhocp1 / rrhocp;
  kappa1 = 1.0;
  kappa2 = kappa1 / rkappa;

  for (current_m = m_val; current_m < n_val; current_m++) {
    
    Delta1 = L0 / (1 << current_m); 

    for (LEVEL = current_m + 1; LEVEL <= n_val; LEVEL++) {
      
      char cmd[256];
      sprintf(cmd, "mkdir -p output/deltax=%d/LEVEL_%d", (1 << current_m), LEVEL);
      system(cmd);

      init_grid (1 << LEVEL);
      TOLERANCE = 1e-4;

      printf("\n--- Lancement: m = %d (Delta1 = 1/%d), LEVEL = %d ---\n", 
             current_m, (1 << current_m), LEVEL);
      
      run();
    }
  }
  return 0;
}

event init (t = 0) {
  char filename[100];
  sprintf (filename, "output/deltax=%d/LEVEL_%d/data.vof", (1 << current_m), LEVEL);
  fpout = fopen (filename, "w");

  foreach() {
    double r1 = sqrt(sq(x) + sq(y));
    d[] =  r1 - R0;
  }

  foreach() {
    if (d[] >= Delta1/2.)
      sf[] = 0.;
    else if (d[] <= -Delta1/2.)
      sf[] = 1.;
    else
      sf[] = 0.5 - d[]/Delta1;
  }

  foreach() {
    fh1[] = rhocp1 * sf[];
    fh2[] = 0.0;
  }
  
}

event stability (i++) {
  double dx = L0 / (1 << LEVEL);
  double kappamax = max(kappa1, kappa2);
  double dtdiff = 0.25 * sq(dx)/kappamax;
  dt = dtnext (dtdiff);
}

event diffusion (i++) {

  scalar T[], s[], lambda[];
  foreach() {
    double rhocp = rhocp1*sf[] + rhocp2*(1. - sf[]);
    T[] = (fh1[] + fh2[])/rhocp;
    s[] = -rhocp*T[];
    lambda[] = -rhocp;
  }

  face vector kappa[];
  foreach_face() {
    double ff = (sf[1] + sf[])/2.;
    kappa.x[] = (kappa1*ff + kappa2*(1. - ff))*dt;
  }

  TOLERANCE = 1e-4;
  poisson (T, s, kappa, lambda);

  foreach() {
    fh1[] = rhocp1*T[]*sf[];
    fh2[] = rhocp2*T[]*(1. - sf[]);
  }
}


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({sf}, (double[]){1.e-2}, LEVEL);
}
#endif

scalar T[], T1[], T2[];
event output (i++) {
  foreach() {
    double rhocp = rhocp1*sf[] + rhocp2*(1. - sf[]);
    T[] = (fh1[] + fh2[])/rhocp;
  }
  stats sT = statsf (T);

  foreach() {
    T1[] = fh1[]/rhocp1*sf[];
    T2[] = fh2[]/rhocp2*(1. - sf[]);
  }

  printf ("%g %g %g \n", t, sT.max, sT.min);
  fprintf (fpout, "%g %g %g \n", t, sT.max, sT.min);

  if (sT.max - sT.min < 0.01) {
    if (fpout) {
      fclose(fpout);
      fpout = NULL;
    }
    return 1; 
  }
}

FILE * fp;
event video (i += 50; t <= tend) {

  char sf_mp4[100], T_mp4[100];
  sprintf (sf_mp4, "output/deltax=%d/LEVEL_%d/sf.mp4", (1 << current_m), LEVEL);
  sprintf (T_mp4,  "output/deltax=%d/LEVEL_%d/T.mp4", (1 << current_m), LEVEL);

  output_ppm (sf, fp, n=2048, min=0, max=1, file=sf_mp4);
  output_ppm (T,  fp, n=2048, min=0, max=1, file=T_mp4);
}

event end_sim (t = tend) {
  if (fpout) {
    fclose(fpout);
    fpout = NULL;
  }
}