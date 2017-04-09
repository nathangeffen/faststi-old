#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main()
{
  const gsl_rng_type * T;
  gsl_rng * r;

  double mu = 9.66430287649229;
  double sigma = 1.25024470140129;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (int i = 0; i < 100; i++)  {
    double x = gsl_ran_weibull (r, mu, sigma);
    printf ("%.2f\n", x);
  }

  gsl_rng_free (r);
  return 0;
}
