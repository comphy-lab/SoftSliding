#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];
double vcm1 , vcm2;

scalar f1[];
scalar f2[];
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);

  // boundary conditions
  u.t[left] = dirichlet(0.0);
  f1[left] = dirichlet(0.0);
  f2[left] = dirichlet(1.0);
  restore (file = filename);
  f1.prolongation = fraction_refine;
  boundary((scalar *){f1, u.x, u.y});

  double sumv1 = 0.;
  double sumv2 = 0.;
  double sumf = 0.;

  foreach() {
    sumv1 += clamp(f1[], 0., 1.)*u.x[];
    sumv2 += clamp(f1[], 0., 1.)*u.y[];

    sumf += clamp(f1[], 0., 1.);    
  }
  
  vcm1 = sumv1/sumf;
  vcm2 = sumv2/sumf;

  boundary((scalar *){f1, u.x, u.y});

  FILE * fp = ferr;
  fprintf(fp, "%f %f %f\n", vcm1, vcm2, t);

  fflush (fp);
  fclose (fp);
}