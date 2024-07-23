#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f1[];
char filename[80];

int main(int a, char const *arguments[]){
  f1[left] = dirichlet(0.);
  // printf("hello");
  sprintf(filename, "%s", arguments[1]);
  // printf("Filename: %s\n", filename);

  restore (file = filename);
  f1.prolongation = fraction_refine;
  boundary((scalar *){f1});

  FILE * fp = ferr;
  output_facets(f1,fp);
  fflush (fp);
  fclose (fp);
}