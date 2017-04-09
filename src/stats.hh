#ifndef __SIM_STATS_H__
#define __SIM_STATS_H__

#include <random>
#include <iostream>

namespace sim  {

  template <typename T>
  T mean(const std::vector<T> & values)
  {
    T total = 0;
    for (auto & v : values)
      total += v;
    return total / (T) values.size();
  }

  template <typename T>
  T median(std::vector<T> & values)
  {
    T answer;
    size_t middle;

    if (values.size() == 0)
      return 0.0;

    if (values.size() == 1)
      return values[0];

    middle = values.size() / 2;
    nth_element(values.begin(), values.begin() + middle, values.end());

    if (values.size() % 2 == 0) {
      answer = *(values.begin() + middle) +
	*max_element(values.begin(), values.begin() + middle);
      answer /= 2.0;
    } else {
      answer = *(values.begin() + middle);
    }
    return answer;
  }


  template <typename RealType = double>
  class beta_distribution
  {
  public:
    typedef RealType result_type;

    class param_type
    {
    public:
      typedef beta_distribution distribution_type;

      explicit param_type(RealType a = 2.0, RealType b = 2.0)
	: a_param(a), b_param(b) { }

      RealType a() const { return a_param; }
      RealType b() const { return b_param; }

      bool operator==(const param_type& other) const
      {
	return (a_param == other.a_param &&
		b_param == other.b_param);
      }

      bool operator!=(const param_type& other) const
      {
	return !(*this == other);
      }

    private:
      RealType a_param, b_param;
    };

    explicit beta_distribution(RealType a = 2.0, RealType b = 2.0)
      : a_gamma(a), b_gamma(b) { }
    explicit beta_distribution(const param_type& param)
      : a_gamma(param.a()), b_gamma(param.b()) { }

    void reset() { }

    param_type param() const
    {
      return param_type(a(), b());
    }

    void param(const param_type& param)
    {
      a_gamma = gamma_dist_type(param.a());
      b_gamma = gamma_dist_type(param.b());
    }

    template <typename URNG>
    result_type operator()(URNG& engine)
    {
      return generate(engine, a_gamma, b_gamma);
    }

    template <typename URNG>
    result_type operator()(URNG& engine, const param_type& param)
    {
      gamma_dist_type a_param_gamma(param.a()),
	b_param_gamma(param.b());
      return generate(engine, a_param_gamma, b_param_gamma);
    }

    template <typename URNG>
    result_type operator()(double from, double to, URNG& engine)
    {
      double range = to - from;
      return range * generate(engine, a_gamma, b_gamma) + from;
    }


    result_type min() const { return 0.0; }
    result_type max() const { return 1.0; }

    result_type a() const { return a_gamma.alpha(); }
    result_type b() const { return b_gamma.alpha(); }

    bool operator==(const beta_distribution<result_type>& other) const
    {
      return (param() == other.param() &&
	      a_gamma == other.a_gamma &&
	      b_gamma == other.b_gamma);
    }

    bool operator!=(const beta_distribution<result_type>& other) const
    {
      return !(*this == other);
    }

  private:
    typedef std::gamma_distribution<result_type> gamma_dist_type;

    gamma_dist_type a_gamma, b_gamma;

    template <typename URNG>
    result_type generate(URNG& engine,
			 gamma_dist_type& x_gamma,
			 gamma_dist_type& y_gamma)
    {
      result_type x = x_gamma(engine);
      return x / (x + y_gamma(engine));
    }
  };
}

#endif
