////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2015 Clemson University.
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
//! \file hydro_helper.cpp
//! \brief A file to define useful hydro helper routines.
//!
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include "hydro_helper.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

typedef std::vector< double > state_type;

//##############################################################################
// Functions.
//##############################################################################

double
acceleration( nnt::Zone& zone, const state_type& x, const double time )
{

  return
    x[1] / ( 3. * zone.getProperty<double>( nnt::s_TAU ) );

}

double rho_function( nnt::Zone& zone, const state_type& x )
{

  return
    zone.getProperty<double>( nnt::s_RHO_0 ) / gsl_pow_3( x[0] );

}

double d_ln_rho_dt_function(
  nnt::Zone& zone,
  const state_type& x
)
{

  return
    rho_function( zone, x ) / zone.getProperty<double>( nnt::s_TAU );

}

void
evolve_function(
  nnt::Zone& zone,
  Libnucnet__NetView * p_view,
  const double d_dt
)
{

  if(
    Libnucnet__Reac__getNumberOfReactions(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_view )
      )
    ) != 0 && d_dt > 1.e-50
  )
  {
    user::safe_evolve( zone, d_dt, d_dt );
  }

} 

double
compute_entropy(
  nnt::Zone& zone
)
{

  return
    compute_thermo_quantity(
      zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    );

}
    
double
t9_from_entropy_root(
  double d_x,
  nnt::Zone& zone,
  Libnucnet__NetView * p_view
)
{

  double d_result;

  gsl_vector * p_abundances, * p_abundance_changes;

  zone.updateProperty(
    nnt::s_T9,
    d_x
  );

  p_abundances =
    Libnucnet__Zone__getAbundances(
      zone.getNucnetZone()
    );

  p_abundance_changes =
    Libnucnet__Zone__getAbundanceChanges(
      zone.getNucnetZone()
    );

  evolve_function(
    zone, 
    p_view,
    zone.getProperty<double>( nnt::s_DTIME )
  ); 

  d_result =
    compute_entropy( zone )
    -      
    zone.getProperty<double>( nnt::s_ENTROPY_PER_NUCLEON );

  Libnucnet__Zone__updateAbundances(
    zone.getNucnetZone(),
    p_abundances
  );

  Libnucnet__Zone__updateAbundanceChanges(
    zone.getNucnetZone(),
    p_abundance_changes
  );

  gsl_vector_free( p_abundances );
  gsl_vector_free( p_abundance_changes );

  return d_result;

}
    
double t9_function( nnt::Zone& zone, Libnucnet__NetView * p_view )
{

  double t9 =
    nnt::compute_1d_root(
      boost::bind(
        t9_from_entropy_root,
        _1,
        boost::ref( zone ),
        p_view
      ),
      zone.getProperty<double>( nnt::s_T9 ),
      1.02
    );

  return t9;

}

void
observer_function(
  nnt::Zone& zone,
  const state_type& x,
  const state_type& dxdt,
  const double d_t
)
{

  double d_dt =
    d_t - zone.getProperty<double>( nnt::s_TIME );

  std::cout <<
    boost::format( "t = %.5e dt = %.5e\n" ) %
    d_t %
    d_dt;

  std::cout <<
    boost::format( "x = {%.5e, %.5e, %.5e}\n" ) %
    x[0] %
    x[1] %
    x[2]; 

  std::cout <<
    boost::format( "dxdt = {%.5e, %.5e, %.5e}\n\n" ) %
    dxdt[0] %
    dxdt[1] %
    dxdt[2];

}

}  // namespace user
