#ifndef RNDM_HPP_
#define RNDM_HPP_

#include "well.hpp"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

typedef boost::math::policies::policy<boost::math::policies::overflow_error<boost::math::policies::errno_on_error> > my_policy;

BOOST_MATH_DECLARE_SPECIAL_FUNCTIONS(my_policy)

typedef boost::math::lognormal_distribution<double, my_policy> lognormal;
typedef boost::math::gamma_distribution<double, my_policy> gamma;

using namespace std;

struct Rndm
{
    Rndm(unsigned int seed);
    ~Rndm();

    double Uniform();
    double Uniform(double under, double upper);
    double Normal(double mean, double sd);
    double Quantile_Beta(double alpha, double beta);

    private:
};

#endif
