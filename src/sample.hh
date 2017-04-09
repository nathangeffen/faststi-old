#ifndef SAMPLE_HH
#define SAMPLE_HH

#include <iostream>
#include <random>
#include <vector>
extern "C" {
#include "ransampl.h"
}

using namespace std;

class Sample {
public:
  Sample() {};
  Sample(vector<double> prob, mt19937* rng) {
    init(prob, rng);
  };
  void init(vector<double> prob, mt19937* rng) {
    rng_ = rng;
    ws_ = ransampl_alloc(prob.size());
    if (!ws_) {
      std::cerr << "Problem allocating sample." << std::endl;
      exit(1);
    }
    double* p = &prob[0];
    ransampl_set(ws_, p);
  };
  unsigned operator()()
  {
    return (unsigned) ransampl_draw(ws_, uni_(*rng_), uni_(*rng_));
  };
  ~Sample() {
    if (ws_)
      ransampl_free(ws_);
  };
private:
  std::uniform_real_distribution<double> uni_;
  ransampl_ws* ws_ = NULL;
  mt19937* rng_;
};

#endif
