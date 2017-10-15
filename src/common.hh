#ifndef COMMON_HH
#define COMMON_HH

#include <random>
#include <boost/random/mersenne_twister.hpp>

#define YEAR 1.0
#define DAY (YEAR / 365)
#define MIN_AGE 12
#define MAX_AGE 101
#define AGE_INTERVAL 5
#define NUM_INTERVALS ( MAX_AGE / AGE_INTERVAL + 1 )

typedef boost::random::mt19937 RNG_TYPE;
// typedef std::mt19937 RNG_TYPE;

#define MALE 0
#define FEMALE 1
#define HOMOSEXUAL 0
#define HETEROSEXUAL 1

#define MAX_AGE_CSV 28
#define EPSILON 0.000001

#define RPM 1
#define RKPM 2
#define CSPM 3

#define HEURISTIC_DISTANCE 0
#define TABLE_DISTANCE 1
#define SYPHILIS_DISTANCE 2

#define GRAPH_ACCURACY 10000

extern thread_local RNG_TYPE rng;

/**
   Macro used for unit testing.  Tests if two values are equal, incrementing
   successes if true, else incrementing failures.

   @param x the first value to test
   @param y the second value to test
   @param successes integer to increment if test is successful
   @param failures integer to increment if test fails
*/

#define TESTEQ(x, y, successes, failures)       \
  do {                                          \
    auto _t1 = (x);                             \
    auto _t2 = (y);                             \
    std::string _t3(#x);                        \
    std::string _t4(#y);                        \
    if (_t1 == _t2) {                           \
      cout << "PASS:\t" << _t3 << " == " << _t4 \
           << "\tLine:" << __LINE__ << "\n";    \
      ++successes;                              \
    }                                           \
    else {                                      \
      cout << "FAIL:\t" << _t3 << " == " << _t4 \
           << "\t" << _t1 << " != " << _t2      \
           << "\tLine:" << __LINE__ << "\n";    \
      ++failures;                               \
    }                                           \
  } while(0)

#endif
