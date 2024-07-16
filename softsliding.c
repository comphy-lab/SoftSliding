/* Title: Drop sliding on a soft solid film
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

# Version 0.0
# Updated: Jul 16, 2024
*/

// 1 is drop, 2 is film and 3 is air

#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase-elasticTF.h"
#include "log-conform-elasticTF.h"
#include "tension.h"

// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define VelErr (1e-2) // error tolerances in velocity
#define DissErr (1e-3)  // error tolerances in dissipation                                      
#define AErr (1e-3) // error tolerance in Conformation tensor
#define MINlevel 4 // minimum level
#define tsnap (0.01)

// Initialization!
#define Xdist (1.040)
#define R2Drop(x,y,z) (sq(x - Xdist) + sq(y))
// domain
#define Ldomain 8                                // Dimension of the domain

// boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = dirichlet(0.0);
f2[left] = dirichlet(1.0);

int MAXlevel;
double tmax;
// Drop
double Ohd; 
// Film
double Ohf, hf, Ec; // De is \infty
double Bond, alphaAngle;
// air
#define RhoA (1e-3)
#define OhA (1e-5)

int main(int argc, char const *argv[]) {

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  MAXlevel = 8;
  tmax = 3.0;

  // Drop
  Ohd = 1e0; // <0.000816/sqrt(816*0.017*0.00075) = 0.008>
  // Film
  Ohf = 1e0;
  hf = 1e0;
  Ec = 0.1;
  Bond = 1e0;
  alphaAngle = pi/6.0;

  fprintf(ferr, "Level %d tmax %g. Ohd %g, Ohf %3.2e, hf %3.2f, Ec %3.2f, De infty \n", MAXlevel, tmax, Ohd, Ohf, hf, Ec);

  L0=Ldomain;
  X0=-hf; Y0=-1.5;
  init_grid (1 << (MINlevel));

  // drop
  rho1 = 1.0; mu1 = Ohd; 
  // film
  rho2 = 1.0; mu2 = Ohf; G2 = Ec;
  // air
  rho3 = RhoA; mu3 = OhA;

  f1.sigma = 1.0; f2.sigma = 1.0;

  run();

}

event acceleration(i++){
  face vector av = a;
  foreach_face(x){
    av.x[] -= Bond*cos(alphaAngle);
  }
  foreach_face(y){
    av.y[] += Bond*sin(alphaAngle);
  }
}

event init(t = 0){
  if(!restore (file = "restart")){
    refine((R2Drop(x,y,z) < 1.44) && (level < MAXlevel)); 
    refine((x < 2*Delta) && (x > -2.*Delta) && (level < MAXlevel));
    fraction (f1, 1. - R2Drop(x,y,z));
    fraction (f2, -x);
  }
}

event adapt(i++){
  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);

  adapt_wavelet ((scalar *){f1, f2, KAPPA1, KAPPA2, u.x, u.y, conform_p.x.x, conform_p.x.y, conform_p.y.y},
  (double[]){fErr, fErr, KErr, KErr, VelErr, VelErr, AErr, AErr, AErr}, 
  MAXlevel, MINlevel);
}

// Outputs
event writingFiles (t = 0, t += tsnap; t <= tmax) {
  dump (file = "restart");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0., vcm = 0., wt = 0.;
  foreach (reduction(+:ke), reduction(+:vcm), reduction(+:wt)){
    ke += (0.5*rho(f1[], f2[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    vcm += (f1[]*u.x[])*sq(Delta);
    wt += f1[]*sq(Delta);
  }
  if (wt > 0.0) vcm /= wt;
  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke vcm\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. Ohd %g, Ohf %3.2e, hf %3.2f, Ec %3.2f, De infty \n", MAXlevel, tmax, Ohd, Ohf, hf, Ec);
      fprintf (fp, "i dt t ke vcm\n");
    } else {
      fp = fopen ("log", "a");
    }
    fprintf (fp, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
    fclose(fp);
    fprintf (ferr, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
  }
  assert(ke > -1e-10);
  assert(ke < 1e2);
  // dump(file = "dumpTest");
}
