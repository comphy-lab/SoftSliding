#include "navier-stokes/centered.h"
#include "fractions.h"
#include "curvature.h"

scalar f2[];
char filename[80];
double ymin;

int main(int a, char const *arguments[]){
  f2[left] = dirichlet(1.0);
  sprintf(filename, "%s", arguments[1]);

  restore (file = filename);
  f2.prolongation = fraction_refine;
  boundary((scalar *){f2});

  scalar pos[];
  position (f2, pos, {1,0});

  ymin = statsf(pos).min;

  FILE * fp = ferr;
  fprintf (fp, "%g %g\n", t, -ymin);

  fflush (fp);
  fclose (fp);
}