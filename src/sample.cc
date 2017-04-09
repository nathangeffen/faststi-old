#include "sample.hh"

int main()
{
  mt19937 rng(5);
  std::vector<double> p = {0.1, 0.1, 0.2, 0.4, 0.1, 0.05, 0.05};
  Sample sample(p, &rng);

  vector<unsigned> cumul(p.size(), 0);
  std::uniform_real_distribution<double> uni;
  for (int j = 0; j < 10000; ++j) {
    unsigned i = sample();
    cumul[i] += 1;
  }

  int total = 0;
  for (int i = 0; i < cumul.size(); ++i) {
    cout << cumul[i] << " ";
    total += cumul[i];
  }
  cout << endl;

  for (int i = 0; i < cumul.size(); ++i) {
    cout << (double) cumul[i] / total  << " ";
  }
  cout << endl;

  return 0;
}
