#ifndef SAMPLE_HH
#define SAMPLE_HH

/**
   C++ interface to ransampl C code.
 */

#include <iostream>
#include <random>
#include <vector>
extern "C" {
#include "ransampl.h"
}

using namespace std;

template <class RNG=std::mt19937>
class Sample {
public:
  /**
     Default constructor does nothing. Init must then be called to initialize
     the class.
   */
  Sample() {};

  /**
     Constructor that calls init (see init for parameter explanation).
   */
  Sample(vector<double> prob, RNG* rng) {
    init(prob, rng);
  };

  /**
     Initializes the sampling mechanism.

     @param prob[in] Weighted probability vector
     @param rng[in,out] Random number generator
   */
  void init(vector<double> prob, RNG* rng) {
    rng_ = rng;
    ws_ = ransampl_alloc(prob.size());
    if (!ws_) {
      std::cerr << "Problem allocating sample." << std::endl;
      exit(1);
    }
    double* p = &prob[0];
    ransampl_set(ws_, p);
  };

  /**
     Draws a random value from weighted probability vector.
   */
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
  RNG* rng_;
};

#endif
