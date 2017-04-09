#include <cstdio>
#include <random>

int main()
{
  double mu;
  double sigma;
  printf("Enter mu: ");
  scanf("%lf", &mu);
  printf("Enter sigma: ");
  scanf("%lf", &sigma);
  printf("mu: %f sigma: %f\n", mu, sigma);
  std::mt19937 rng;
  std::weibull_distribution<double> dist(sigma, mu);

  for (int i = 0; i < 100; i++)  {
    double x = dist(rng);
    printf ("%.2f\n", x);
  }

  return 0;
}
