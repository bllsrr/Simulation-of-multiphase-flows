void initial_condition (int Nvort, vector u) {


  scalar omega[];
	
  srand(123); 
  int Nimages = 1;
  double G1[Nvort/2];
  double G2[Nvort/2];
  double avort = 10.0/pow(2,7); //third run
  for (int i=0; i < Nvort/2; i++) {
    G1[i] = noise();
    G2[i] = noise();
  }

  double Gsum1 = 0.;
  double Gsum2 = 0.;
  for (int i=0; i < Nvort/2; i++) {
	Gsum1 += G1[i];
	Gsum2 += G2[i];
  }
  double acorr = Gsum2/(Gsum1 + Gsum2);
  

  for (int i=0; i < Nvort/2; i++) {
    double x0 = 0.5*noise();
    double y0 = 0.25 + (0.25)*noise();
    if ( fabs(y0) - 0.25 < avort)
      y0 += avort*pow(-1,i);

    double amp0 = G1[i]*acorr;
    G1[i] = G1[i]*acorr*M_PI*pow(avort,2);
    foreach() {
      for (int j=-Nimages; j <= Nimages; j++) 
        for (int k=-Nimages; k <= Nimages; k++) 
          omega[] += amp0*exp( -( pow(x-x0-k*L0,2) + pow(y-y0-j*L0,2))/pow(avort,2) );
    }
  }

  for (int i=0; i < Nvort/2; i++) {
    double x0 = 0.5*noise();
    double y0 = -0.25 + (0.25)*noise();
    if ( fabs(y0) - 0.25 < avort)
      y0 += avort*pow(-1,i);

    double amp0 = G2[i]*(acorr-1.);
    G2[i] = G2[i]*(acorr - 1.)*M_PI*pow(avort,2);
    foreach() {
      for (int j=-Nimages; j <= Nimages; j++) 
        for (int k=-Nimages; k <= Nimages; k++) 
          omega[] += amp0*exp( -( pow(x-x0-k*L0,2) + pow(y-y0-j*L0,2))/pow(avort,2) );
    }
  }

  scalar psi[];
  psi[right]  = dirichlet(0);
  psi[left]   = dirichlet(0);
  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);

  poisson (psi, omega);

  double unorm = 0.;
  
  struct { double x, y; } ff = {-1.,1.};
  foreach() 
    foreach_dimension () {
      u.x[] = ff.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
      unorm += sq(u.x[])*sq(Delta);
    }
  
   foreach() 
    foreach_dimension () 
      u.x[] /= sqrt(unorm/sq(L0));
  
}
