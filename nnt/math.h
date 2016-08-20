////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2015 Clemson University.
//
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for a variety of useful math utilities.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_MATH_H
#define NNT_MATH_H

#include <iostream>

#include <WnMatrix.h>
#include <Libnucnet.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <boost/math/tools/roots.hpp>
#include <boost/tuple/tuple.hpp>

#include "nnt/param_defs.h"
#include "nnt/wrappers.hpp"

namespace nnt
{

//##############################################################################
// Prototypes.
//##############################################################################

double
compute_1d_root(
  boost::function<double(double)>,
  double,
  double,
  double
);
  
double
compute_1d_root(
  boost::function<double(double)>,
  double,
  double
);

std::pair<double,double>
default_guess_function(
  nnt::Zone&,
  const std::string&,
  double
);
  
std::pair<double,double>
default_log10_guess_function(
  nnt::Zone&,
  const std::string&,
  double
);
  
std::pair<double,double>
bilinear_interpolation(
  gsl_vector *, gsl_vector *, gsl_matrix *, double, double
);

std::pair<double,double>
two_d_interpolation(
  gsl_vector *, gsl_vector *, gsl_matrix *, double, double
);

double
linear_interpolation(
  gsl_vector *, gsl_vector *, double
);

int
get_table_position(
  gsl_vector *, gsl_vector *, double, double
);

double
compute_derivative(
  boost::function<double(double)>,
  const double&,
  double
);

double
compute_derivative(
  boost::function<double(double)>,
  const double&
);

double
linear_interpolation( gsl_vector *, gsl_vector *, double );

double
spline_interpolation( gsl_vector *, gsl_vector *, double );

} // namespace nnt

#endif // NNT_MATH_H
