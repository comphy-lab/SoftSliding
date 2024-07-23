#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f2[];
char filename[80];

int main(int a, char const *arguments[]){
  f2[left] = dirichlet(1.0);
  // printf("hello");
  sprintf(filename, "%s", arguments[1]);
  // printf("Filename: %s\n", filename);

  restore (file = filename);
  f2.prolongation = fraction_refine;
  boundary((scalar *){f2});

  FILE * fp = ferr;
  output_facets(f2,fp);
  fflush (fp);
  fclose (fp);
}