//------------------------------------------------------------------------------
// BOOST used for statistical distributions NOT USED October 2022
// return NAN for invalid inputs
// see https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/stat_tut/weg/error_eg.html
// and https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/changing_policy_defaults.html
#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
// must follow define domain error policy...
#include <boost/math/distributions.hpp>     
//------------------------------------------------------------------------------

// RcppArmadillo.h pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
