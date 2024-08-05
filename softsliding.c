/* Title: Drop sliding on a soft solid film
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

# Version 1.0
# Updated: Aug 5, 2024

# changelog Aug 5, 2024
* This code uses the reduced gravity formulation described here: (https://www.annualreviews.org/content/journals/10.1146/annurev-fluid-122316-045034)[https://www.annualreviews.org/content/journals/10.1146/annurev-fluid-122316-045034]
*/

// 1 is drop, 2 is film and 3 is air

#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase-nonCoalescing-elastic.h"
#include "log-conform-elastic.h"
#include "tension.h"
#include "reduced-three-phase-nonCoalescing.h"

// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define VelErr (1e-2) // error tolerances in velocity
#define AErr (1e-3) // error tolerance in Conformation tensor
#define MINlevel 4 // minimum level
#define tsnap (0.01)

// Initialization!
// this will be done using the restart file only!
// boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = dirichlet(0.0);
f2[left] = dirichlet(1.0);

u.n[right] = neumann(0.0);
p[right] = dirichlet(0.0);

int MAXlevel;
double tmax;
// Drop
double Ohd; 
// Film
double Ohf, hf, Ec; // De is \infty
double Bond, alphaAngle;
double Ldomain;
// air
#define RhoA (1e-3)
#define OhA (1e-5)

int main(int argc, char const *argv[]) {

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  MAXlevel = atoi(argv[1]);
  tmax = atof(argv[2]);

  // Drop
  Ohd = 1e0; // <0.000816/sqrt(816*0.017*0.00075) = 0.008>
  // Film
  Ohf = 1e0;
  hf = atof(argv[3]); // ratio of the film thickness to the drop radius, log scale: 0.01--1 or so
  Ec = atof(argv[4]); // Elasto-capillary number: 1e-4 (very soft) to 1e3 (very stiff)
  Bond = 1e0; // Bond number: we will keep this fixed
  alphaAngle = pi*atof(argv[5])/180; // inclination angle of the drop: user should define in degrees. 10-60 degrees for the initial runs. 
  Ldomain = 32.0; // Dimension of the domain: should be large enough to get a steady solution to drop velocity.

  fprintf(ferr, "Level %d tmax %g. Ohd %g, Ohf %3.2e, hf %3.2f, Ec %3.2f, Bo %3.2f, alpha %3.2f, De infty \n", MAXlevel, tmax, Ohd, Ohf, hf, Ec, Bond, alphaAngle);

  L0=Ldomain;
  X0=-hf; Y0=-2.0;
  init_grid (1 << (9));

  // drop
  rho1 = 1.0; mu1 = Ohd; G1 = 0.;
  // film
  rho2 = 1.0; mu2 = Ohf; G2 = Ec;
  // air
  rho3 = RhoA; mu3 = OhA; G3 = 0.;

  f1.sigma = 1.0; f2.sigma = 1.0;

  // for sliding
  Bf1.x = -Bond*cos(alphaAngle);
  Bf2.x = -Bond*cos(alphaAngle);

  Bf1.y = -Bond*sin(alphaAngle);
  Bf2.y = -Bond*sin(alphaAngle);

  run();

}


event init(t = 0){
  if (!restore (file = "restart")) {
    fprintf(ferr, "Restart file not found. Trying resetting the system to its base state.\n");
    if(!restore (file = "equilibriumSolution")){
      fprintf(ferr, "Equilibrium solution not found. You should run the equilibrium code first. Exiting.\n");
      return 1;
    }
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
    vcm += (f1[]*u.y[])*sq(Delta);
    wt += f1[]*sq(Delta);
  }
  if (wt > 0.0) vcm /= wt;
  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke vcm\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. Ohd %g, Ohf %3.2e, hf %3.2f, Ec %3.2f, Bo %3.2f, alpha %3.2f, De infty \n", MAXlevel, tmax, Ohd, Ohf, hf, Ec, Bond, alphaAngle);
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
