////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
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
//
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A file containing useful math routines
//!
////////////////////////////////////////////////////////////////////////////////

#include "nnt/math.h"

namespace nnt
{

//############################################################################
// bilinear_interpolation()
//############################################################################

std::pair<double,double>
bilinear_interpolation(
  gsl_vector *p_x1,
  gsl_vector *p_x2,
  gsl_matrix *p_matrix,
  double d_x1,
  double d_x2
)
{

  size_t j, k = 0;
  double d_y1 = 0, d_y2 = 0, d_y3 = 0, d_y4 = 0, d_t = 0, d_u = 0;
  double d_diff;

  if(
   d_x1 < gsl_vector_get( p_x1, 0L ) ||
   d_x1 > gsl_vector_get( p_x1, WnMatrix__get_gsl_vector_size( p_x1 ) - 1 )
  )
  {
    fprintf( stderr, "x1 out of range.\n" );
    exit( EXIT_FAILURE );
  }

  if(
   d_x2 < gsl_vector_get( p_x2, 0L ) ||
   d_x2 > gsl_vector_get( p_x2, WnMatrix__get_gsl_vector_size( p_x2 ) - 1 )
  )
  {
    fprintf( stderr, "x2 out of range.\n" );
    exit( EXIT_FAILURE );
  }

  j =
    gsl_interp_bsearch(
      WnMatrix__get_gsl_vector_array( p_x1 ),
      d_x1,
      0L,
      WnMatrix__get_gsl_vector_size( p_x1 ) - 1
    );

  k =
    gsl_interp_bsearch(
      WnMatrix__get_gsl_vector_array( p_x2 ),
      d_x2,
      0L,
      WnMatrix__get_gsl_vector_size( p_x2 ) - 1
    );

  d_y1 = gsl_matrix_get( p_matrix, j, k );
  d_y2 = gsl_matrix_get( p_matrix, j+1, k );
  d_y3 = gsl_matrix_get( p_matrix, j+1, k+1 );
  d_y4 = gsl_matrix_get( p_matrix, j, k+1 );

  d_diff = fabs( d_y1 - d_y2 ); 
  d_diff = GSL_MAX( d_diff, fabs( d_y1 - d_y3 ) );
  d_diff = GSL_MAX( d_diff, fabs( d_y1 - d_y4 ) );
  d_diff = GSL_MAX( d_diff, fabs( d_y2 - d_y3 ) );
  d_diff = GSL_MAX( d_diff, fabs( d_y2 - d_y4 ) );
  d_diff = GSL_MAX( d_diff, fabs( d_y3 - d_y4 ) );

  d_t =
    ( d_x1 - gsl_vector_get( p_x1, j ) ) /
    ( gsl_vector_get( p_x1, j+1 ) - gsl_vector_get( p_x1, j ) );
  d_u =
    ( d_x2 - gsl_vector_get( p_x2, k ) ) /
    ( gsl_vector_get( p_x2, k+1 ) - gsl_vector_get( p_x2, k ) );

  return
    std::make_pair(
      (1. - d_t ) * (1. - d_u ) * d_y1 +
      d_t * (1. - d_u ) * d_y2 +
      d_t * d_u * d_y3 +
      (1. - d_t ) * d_u * d_y4,
      d_diff
    );

}

//############################################################################
// linear_interpolation()
//############################################################################

/**
  \brief Routine to do a linear interpolation.

  \param p_x A pointer to a gsl_vector giving the array of independent
             variables.
  \param p_y A pointer to a gsl_vector giving the array of dependent variables.
  \param d_x The independent variable value at which to do the interpolation.
  \return The interpolated value.  If d_x is smaller than the smallest
          value in the input array p_x, the value of p_y corresponding to
          the smallest value of p_x is returned.
          If d_x is larger than the largest value in p_x, the value of p_y
          corresponding to the largest value of p_x is returned.
*/

double
linear_interpolation( gsl_vector *p_x, gsl_vector *p_y, double d_x )
{

  size_t j;

  if( d_x < WnMatrix__get_gsl_vector_array( p_x )[0] )
    return
      gsl_vector_get( p_y, 0L );
  else if(
    d_x >=
    gsl_vector_get(
      p_x,
      WnMatrix__get_gsl_vector_size( p_x ) - 1
    )
  )
    return
      gsl_vector_get(
	p_y,
	WnMatrix__get_gsl_vector_size( p_x ) - 1
      );
  else
  {
    j =
      gsl_interp_bsearch(
	WnMatrix__get_gsl_vector_array( p_x ),
	d_x,
	0L,
	WnMatrix__get_gsl_vector_size( p_x ) - 1
      );
    return
      gsl_vector_get( p_y, j ) +
      ( gsl_vector_get( p_y, j + 1 ) - gsl_vector_get( p_y, j ) ) *
      ( d_x - gsl_vector_get( p_x, j ) ) /
      ( gsl_vector_get( p_x, j + 1 ) - gsl_vector_get( p_x, j ) );
  }

} 

//############################################################################
// two_d_interpolation()
//############################################################################

/**
  \brief Routine to do a two_d interpolation.

  \param p_x1 A pointer to a gsl_vector giving the array of the first
             independent variables.
  \param p_x2 A pointer to a gsl_vector giving the array of the second
             independent variables.
  \param p_matrix A pointer to a gsl_matrix giving the dependent variables
             corresponding to the independent variabls in p_x1 and p_x2.
  \param d_x1 The first independent variable value at which to do the
             interpolation.
  \param d_x2 The second independent variable value at which to do the
             interpolation.
  \return A pair.  The first value is the interpolated value while
           the second is the maximum difference between interpolation
           points.  If d_x1 and/or d_x2 lie outside the table,
           the routine extrapolates from the available data.  If the
           data lie outside the table, the second of the pair is returned
           as 9999.
*/

std::pair<double,double>
two_d_interpolation(
  gsl_vector *p_x1,
  gsl_vector *p_x2,
  gsl_matrix *p_matrix,
  double d_x1,
  double d_x2
)
{

  gsl_vector_view view;
  int i_table_position;

  i_table_position = get_table_position( p_x1, p_x2, d_x1, d_x2 );

  switch( i_table_position )
  {

    case 0:
      return
        std::make_pair( gsl_matrix_get( p_matrix, 0L, 0L ), 9999. );

    case 1:
      return 
        std::make_pair(
	  gsl_matrix_get(
	    p_matrix,
	    0L,
	    WnMatrix__get_gsl_vector_size( p_x2 ) - 1
	  ),
          9999.
        ); 

    case 2:
      view = gsl_matrix_row( p_matrix, 0L );
      return
        std::make_pair(
          linear_interpolation( p_x2, &view.vector, d_x2 ),
          9999.
        );

    case 3:
      return
        std::make_pair(
  	  gsl_matrix_get(
	    p_matrix,
	    WnMatrix__get_gsl_vector_size( p_x1 ) - 1,
	    0L
	  ),
          9999.
        );

    case 4:
      return
        std::make_pair(
          gsl_matrix_get(
	    p_matrix,
	    WnMatrix__get_gsl_vector_size( p_x1 ) - 1,
	    WnMatrix__get_gsl_vector_size( p_x2 ) - 1
	  ),
          9999.
        ); 

    case 5:
      view =
	gsl_matrix_row(
	  p_matrix,
	  WnMatrix__get_gsl_vector_size( p_x1 ) - 1
	);
      return
        std::make_pair(
 	  linear_interpolation( p_x2, &view.vector, d_x2 ),
          9999.
        );

    case 6:
      view = gsl_matrix_column( p_matrix, 0L );
      return
	std::make_pair(
          linear_interpolation( p_x1, &view.vector, d_x1 ),
          9999.
        );

    case 7:
      view =
	gsl_matrix_column(
	  p_matrix,
	  WnMatrix__get_gsl_vector_size( p_x2 ) - 1
	);
      return
        std::make_pair(
	  linear_interpolation( p_x1, &view.vector, d_x1 ),
          9999.
        );

    case 8:
      return
	bilinear_interpolation(
	  p_x1,
	  p_x2,
	  p_matrix,
	  d_x1,
	  d_x2
	);

    default:
      fprintf( stderr, "Invalid case.\n" );
      exit( EXIT_FAILURE );

  }

}

//############################################################################
// get_table_position()
//############################################################################

int
get_table_position(
  gsl_vector *p_x1,
  gsl_vector *p_x2,
  double d_x1,
  double d_x2
)
{

  if( d_x1 < gsl_vector_get( p_x1, 0L ) )
  {
    if( d_x2 < gsl_vector_get( p_x2, 0L ) )
      return 0;
    else if(
      d_x2 >
      gsl_vector_get(
	p_x2,
	WnMatrix__get_gsl_vector_size( p_x2 ) - 1
      )
    )
      return 1;
    else
      return 2;
  }
  else if(
    d_x1 >
    gsl_vector_get(
      p_x1,
      WnMatrix__get_gsl_vector_size( p_x1 ) - 1
    )
  )
  {
    if( d_x2 < gsl_vector_get( p_x2, 0L ) )
      return 3;
    else if(
      d_x2 >
      gsl_vector_get(
	p_x2,
	WnMatrix__get_gsl_vector_size( p_x2 ) - 1
      )
    )
      return 4;
    else
      return 5;
  }
  else
  {
    if( d_x2 < gsl_vector_get( p_x2, 0L ) )
      return 6;
    else if(
      d_x2 >
      gsl_vector_get(
	p_x2,
	WnMatrix__get_gsl_vector_size( p_x2 ) - 1
      )
    )
      return 7;
    else
      return 8;
  }

}

//##############################################################################
// compute_1d_root().
//##############################################################################

/**
  \brief Routine to find the root of a function.

  \param quantity_func  The quantity function.
  \param guess  The initial guess.
  \param factor The bracketing expansion factor.
  \param d_precision The fraction of the maximum possible precision for the root (optional--the default is 1.0).
  \return The root of the function.
*/

double
compute_1d_root(
  boost::function<double(double)> quantity_func,
  double guess,
  double factor,
  double d_precision
)
{
  
  const boost::uintmax_t max_iter = 100;
  boost::uintmax_t iter;
  bool is_rising;

  if( quantity_func( guess * factor ) > quantity_func( guess / factor ) )
    is_rising = true;
  else
    is_rising = false;

  iter = max_iter;

  int digits = std::numeric_limits<double>::digits;
  int get_digits = static_cast<int>( digits * d_precision );

  boost::math::tools::eps_tolerance<double> tol( get_digits );

  std::pair<double, double> r =
    boost::math::tools::bracket_and_solve_root(
      quantity_func,
      guess,
      factor,
      is_rising,
      tol,
      iter
    );

  if( iter >= max_iter )
  {
    std::cerr << "Root not found--too many iterations." << std::endl;
    exit( EXIT_FAILURE );
  }
    
  return (r.first + r.second) / 2.;

}

double
compute_1d_root(
  boost::function<double(double)> quantity_func,
  double guess,
  double factor
)
{

  return compute_1d_root( quantity_func, guess, factor, 1.0 );

}
 
//##############################################################################
// default_guess_function
//##############################################################################

std::pair<double,double>
default_guess_function(
  nnt::Zone& zone,
  const std::string& s_quantity_key,
  double d_factor
)
{

  double guess;

  if( zone.hasProperty( s_quantity_key ) )
    guess = zone.getProperty<double>( s_quantity_key );
  else
    guess = 1;

  return std::make_pair( guess, d_factor );

}  

//##############################################################################
// default_log10_guess_function
//##############################################################################

std::pair<double,double>
default_log10_guess_function(
  nnt::Zone& zone,
  const std::string& s_quantity_key,
  double d_factor
)
{

  double guess;

  if( zone.hasProperty( s_quantity_key ) )
    guess = log10( zone.getProperty<double>( s_quantity_key ) );
  else
    guess = 1;

  return std::make_pair( guess, d_factor );

}  

//##############################################################################
// compute_derivative().
//##############################################################################

/**
  \brief Routine to find the derivative of a function.

  \param func The function.
  \param x The value at which to compute the derivative.
  \param d_step_fraction The fractional step size to use for the derivative
         (optional--default = 0.01).
  \return The derivative at the input value.
*/

double
compute_derivative(
  boost::function<double(double)> func,
  const double& x,
  double d_step_fraction
)
{

  const double dx = x * d_step_fraction;

  const double m1 = ( func(x + dx) - func(x - dx) ) / 2;
  const double m2 = ( func(x + 2 * dx) - func(x - 2 * dx) ) / 4;
  const double m3 = ( func(x + 3 * dx) - func(x - 3 * dx) ) / 6;

  return ( ( 15 * m1 - 6 * m2 ) + m3) / ( 10 * dx );
  
}

double
compute_derivative(
  boost::function<double(double)> func,
  const double& x
)
{

  return compute_derivative( func, x, 0.01 );

}

//##############################################################################
// spline_interpolation().
//##############################################################################

/**
  \brief Routine to do a spline interpolation.

  \param p_x A pointer to a gsl_vector giving the array of independent
             variables.
  \param p_y A pointer to a gsl_vector giving the array of dependent variables.
  \param d_x The independent variable value at which to do the interpolation.
  \return The interpolated value.  If d_x is smaller than the smallest
          value in the input array p_x, the value of p_y corresponding to
          the smallest value of p_x is returned.
          If d_x is larger than the largest value in p_x, the value of p_y
          corresponding to the largest value of p_x is returned.
*/

double
spline_interpolation( gsl_vector * p_x, gsl_vector * p_y, double d_x )
{

  gsl_interp_accel *p_acc;
  gsl_spline *p_spline;
  double d_result;

  if( d_x < WnMatrix__get_gsl_vector_array( p_x )[0] )
  {
    return
      gsl_vector_get( p_y, 0L );
  }
  else if(
    d_x >=
    gsl_vector_get(
      p_x,
      WnMatrix__get_gsl_vector_size( p_x ) - 1
    )
  )
  {
    return
      gsl_vector_get(
	p_y,
	WnMatrix__get_gsl_vector_size( p_x ) - 1
      );
  }

  p_acc = gsl_interp_accel_alloc();

  p_spline =
    gsl_spline_alloc(
      gsl_interp_cspline,
      WnMatrix__get_gsl_vector_size( p_x )
    );

  gsl_spline_init(
    p_spline,
    WnMatrix__get_gsl_vector_array( p_x ),
    WnMatrix__get_gsl_vector_array( p_y ),
    WnMatrix__get_gsl_vector_size( p_x )
  );

  d_result = gsl_spline_eval( p_spline, d_x, p_acc );

  //============================================================================
  // Clean up and return.
  //============================================================================

  gsl_spline_free( p_spline );
  gsl_interp_accel_free( p_acc );

  return d_result;

}

} // namespace nnt

