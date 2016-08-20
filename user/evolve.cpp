//////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Code for evolving a network.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes. 
//##############################################################################

#include "user/evolve.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */
namespace user
{

//##############################################################################
// evolve()
//##############################################################################

/**
 * \brief Evolve a Nucnet Tools zone over the currently defined time
 *        step (property s_DTIME) for the zone.
 *
 * \param zone A Nucnet Tools zone.
 * \return Number of iterations if successful, -1 if not successful.
 */

int
evolve( 
  nnt::Zone& zone
) {

  WnMatrix *p_matrix; 
  size_t i_iter;
  gsl_vector *p_y_old, *p_rhs, *p_sol, *p_work;
  double d_dt;
  std::pair<double,double> check;
  
  //==========================================================================
  // Evolve NSE + weak rates, if appropriate.
  //==========================================================================

  if(
    zone.hasProperty( nnt::s_EVOLVE_NSE_PLUS_WEAK_RATES ) &&
    zone.getProperty<std::string>( nnt::s_EVOLVE_NSE_PLUS_WEAK_RATES ) == "yes"
  )
  {
    evolve_nse_plus_weak_rates( zone );
    return 1;
  }

  //============================================================================
  // Get timestep. 
  //============================================================================

  d_dt = zone.getProperty<double>( nnt::s_DTIME );

  //============================================================================
  // Save the old abundances.
  //============================================================================

  p_y_old = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  //============================================================================
  // Newton-Raphson Iterations.
  //============================================================================

  for( i_iter = 1; i_iter <= I_ITMAX; i_iter++ ) {

    //--------------------------------------------------------------------------
    // Get matrix and rhs vector.
    //--------------------------------------------------------------------------

    boost::tie( p_matrix, p_rhs ) = get_evolution_matrix_and_vector( zone );

    //--------------------------------------------------------------------------
    // Add 1/dt to diagonal.
    //--------------------------------------------------------------------------

    WnMatrix__addValueToDiagonals(
      p_matrix,
      1.0 / zone.getProperty<double>( nnt::s_DTIME )
    );

    //--------------------------------------------------------------------------
    // Correct vector for iteration.
    //--------------------------------------------------------------------------

    p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );
    gsl_vector_sub( p_work, p_y_old );
    gsl_vector_scale( p_work, 1. / d_dt );
    gsl_vector_sub( p_rhs, p_work );

    gsl_vector_free( p_work );

    //--------------------------------------------------------------------------
    // Solve matrix equation.
    //--------------------------------------------------------------------------

    p_sol = solve_matrix_for_zone( zone, p_matrix, p_rhs );

    //--------------------------------------------------------------------------
    // Check solution.
    //--------------------------------------------------------------------------

    check = check_matrix_solution( zone, p_sol );

    //--------------------------------------------------------------------------
    // Update abundances.
    //--------------------------------------------------------------------------

    p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

    gsl_vector_add( p_work, p_sol );

    Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_work );

    gsl_vector_free( p_work );

    //--------------------------------------------------------------------------
    // Free matrix, p_rhs, and p_sol. Remember
    // Libnucnet__Zone__computeJacobianMatrix returns a new matrix and
    // Libnucnet__computeFlowVector and WnMatrix__solve return new gsl_vectors
    // each time they are called.
    //--------------------------------------------------------------------------

    WnMatrix__free( p_matrix ); 
    gsl_vector_free( p_rhs );
    gsl_vector_free( p_sol );

    //--------------------------------------------------------------------------
    // Exit iterations if converged.
    //--------------------------------------------------------------------------

    if( zone.hasProperty( nnt::s_NEWTON_RAPHSON_CONVERGE ) )
    {
      if(
        check.first < zone.getProperty<double>( nnt::s_NEWTON_RAPHSON_CONVERGE )
      )
        break;
    }
    else
    {
      if( check.first < D_MIN ) break;
    }

    //--------------------------------------------------------------------------
    // Return with negative value if large negative abundances.
    //--------------------------------------------------------------------------

    if( zone.hasProperty( nnt::s_LARGE_NEG_ABUND_THRESHOLD ) )
    {
      if( !is_nonneg_abunds( zone ) ) return -1;
    }
      
  }

  //==========================================================================
  // Update abundance changes.
  //==========================================================================

  p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  gsl_vector_sub( p_work, p_y_old );

  Libnucnet__Zone__updateAbundanceChanges( zone.getNucnetZone(), p_work );

  gsl_vector_free( p_work );
  
  //==========================================================================
  // Free allocated memory and return.
  //==========================================================================

  gsl_vector_free( p_y_old );

  return (int) i_iter;

}

//##############################################################################
// default_safe_evolve_check_function()
//##############################################################################

bool
default_safe_evolve_check_function( nnt::Zone& zone, const double d_x_min )
{

  double d_xsum =
    1. - Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 );

  return ( fabs( d_xsum ) < d_x_min || is_nonneg_abunds( zone ) );

}

//##############################################################################
// safe_evolve()
//##############################################################################

/**
 * \brief Evolve a Nucnet Tools zone over a timestep [0,d_dt] using
 *        initial step d_dt1.  Check the resulting abundances with the
 *        currently defined SAFE_EVOLVE_CHECK_FUNCTION.  If the abundances
 *        do not satisfy the check function, cut the time step until they do.
 *        Then build the time step back up until the full elapse time is d_dt.
 *
 * \param zone A Nucnet Tools zone.
 * \param d_dt1 The initial step to try.  If not supplied, taken as the
 *              current s_DTIME property for the zone.
 * \param d_dt The time step for full evolution.  If not supplied, taken as the
 *              current s_DTIME property for the zone.
 * \param d_dt_min The minimum time step to cut to.  (Optional--default is
 *                 1.e-20.)
 */

void
safe_evolve(
  nnt::Zone& zone,
  double d_dt1,
  const double d_dt,
  const double d_dt_min
)
{

  gsl_vector *p_y, *p_y_old;
  double d_t;
  boost::function<bool( nnt::Zone& )> check_f;

  Libnucnet__Zone * p_zone = zone.getNucnetZone();

  p_y_old = Libnucnet__Zone__getAbundances( p_zone );

  zone.updateProperty(
    nnt::s_DTIME,
    d_dt1
  );

  evolve( zone );

  if( !zone.hasFunction( nnt::s_SAFE_EVOLVE_CHECK_FUNCTION ) )
  {
    zone.updateFunction(
      nnt::s_SAFE_EVOLVE_CHECK_FUNCTION,
      static_cast<boost::function<bool( nnt::Zone& )> >(
        boost::bind( default_safe_evolve_check_function, _1, D_X_EPS )
      )
    );
  }

  check_f =
    boost::any_cast<boost::function<bool( nnt::Zone& )> >(
      zone.getFunction( nnt::s_SAFE_EVOLVE_CHECK_FUNCTION )
    );

  if( !check_f( zone ) )
  {

    while( !check_f( zone ) && d_dt1 > d_dt_min )
    {
      Libnucnet__Zone__updateAbundances( p_zone, p_y_old );
      d_dt1 /= 10;
      zone.updateProperty(
        nnt::s_DTIME,
        d_dt1
      );
      evolve( zone );
    }

  }

  d_t = d_dt1;

  while( d_t < d_dt )
  {
    evolve( zone );
    d_t += d_dt1;
    d_dt1 *= 1.15;
    if( d_dt1 + d_t > d_dt ) d_dt1 = d_dt - d_t + 1.e-30;
    zone.updateProperty(
      nnt::s_DTIME,
      d_dt1
    );
  }

  p_y = Libnucnet__Zone__getAbundances( p_zone );

  gsl_vector_sub( p_y, p_y_old );

  Libnucnet__Zone__updateAbundanceChanges( p_zone, p_y );

  zone.updateProperty(
    nnt::s_DTIME,
    d_dt
  );

  gsl_vector_free( p_y );
  gsl_vector_free( p_y_old );

}

void
safe_evolve(
  nnt::Zone& zone,
  double d_dt1,
  const double d_dt
)
{
  safe_evolve( zone, d_dt1, d_dt, 1.e-20 );
}

void
safe_evolve( nnt::Zone& zone )
{

  safe_evolve(
    zone,
    zone.getProperty<double>( nnt::s_DTIME ),
    zone.getProperty<double>( nnt::s_DTIME )
  );

}

//##############################################################################
// is_nonneg_abunds().
//##############################################################################

bool
is_nonneg_abunds( nnt::Zone& zone )
{

  double d_abund_min = 0.;

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      )
    );

  if( zone.hasProperty( nnt::s_LARGE_NEG_ABUND_THRESHOLD ) )
  {
    d_abund_min = zone.getProperty<double>( nnt::s_LARGE_NEG_ABUND_THRESHOLD );
  }

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    double d_abund =
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      );

    if( d_abund < 0 )
    {
      if( fabs( d_abund ) < d_abund_min )
      {
        Libnucnet__Zone__updateSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies(),
          0.
        );
      }
      else
      {
        return false;
      }
    }

  }

  return true;

}

//##############################################################################
// evolve_zone().
//##############################################################################

int
evolve_zone( nnt::Zone& zone, double d_t_end )
{

  Libnuceq *p_equil;
  gsl_vector *p_abundances;
  double d_t = 0., d_dt, d_t9, d_rho; 
  int i_steps = 0;

  //============================================================================
  // Initializations.
  //============================================================================

  std::cout <<
    boost::format( "Start Zone %s with T9 = %s\n" ) %
    Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) %
    zone.getProperty<std::string>( nnt::s_T9_0 ).c_str();

  zone.updateProperty(
    nnt::s_TIME,
    d_t
  );

  d_t9 = zone.getProperty<double>( nnt::s_T9_0 );

  d_rho = zone.getProperty<double>( nnt::s_RHO_0 );


  //============================================================================
  // Normalize abundances so that Xsum = 1.
  //============================================================================

  nnt::normalize_zone_abundances( zone );

  //============================================================================
  // Set to NSE if t9 > 7.
  //============================================================================

  if( d_t9 > 7. )
  {

    d_t9 = 7.;

    d_rho =
      gsl_pow_3(
        d_t9 /
        zone.getProperty<double>( nnt::s_T9_0 )
      ) *
      zone.getProperty<double>( nnt::s_RHO_0 );

    d_t =
      zone.getProperty<double>( nnt::s_TAU )
      *
      log( zone.getProperty<double>( nnt::s_RHO_0 ) / d_rho );
        
    p_equil =
      Libnuceq__new(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      );

    Libnuceq__setYe(
      p_equil, 
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
    );

    Libnuceq__computeEquilibrium(
      p_equil, 
      d_t9,
      d_rho
    );

    p_abundances = Libnuceq__getAbundances( p_equil );

    Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abundances );

    gsl_vector_free( p_abundances );

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

    zone.updateProperty(
      nnt::s_DTIME,
      1.e-6
    );


  }

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  while( d_t < d_t_end )
  {

    i_steps++;

    d_dt = zone.getProperty<double>( nnt::s_DTIME );
    d_t += d_dt;

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

    d_t9 =
      zone.getProperty<double>( nnt::s_T9_0 ) *
      exp(
        -d_t /
        (
          3. * zone.getProperty<double>( nnt::s_TAU )
        )
      );

    if( d_t9 < 1.e-6 ) d_t9 = 1.e-6;

    zone.updateProperty(
      nnt::s_T9,
      d_t9
    );

    d_rho = zone.getProperty<double>( nnt::s_RHO_0 ) *
      exp( -d_t / zone.getProperty<double>( nnt::s_TAU ) );

    if( d_rho < 1.e-18 ) d_rho = 1.e-18;

    zone.updateProperty(
      nnt::s_RHO,
      d_rho
    );

  //============================================================================
  // Evolve system.
  //============================================================================

    if( !evolve( zone ) )
    {
      std::cerr << "Problem converging." << std::endl;
      exit( EXIT_FAILURE );
    }

  //============================================================================
  // Update timestep.
  //============================================================================

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      0.15,
      0.15,
      1.e-10
    );

    if ( d_t + d_dt > d_t_end ) {

      d_dt = d_t_end - d_t;

    }

    zone.updateProperty(
      nnt::s_DTIME,
      d_dt
    );

    zone.updateProperty(
      nnt::s_STEPS,
      i_steps
    );

  //============================================================================
  // Limit the network.
  //============================================================================

    limit_evolution_network( zone );

  }  

  std::cout <<
    boost::format( 
      "                          End Zone %s, %d steps, 1 - xsum = %g\n"
    ) %
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) %
      i_steps %
      ( 1. - Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) );

  return 0;

}

//##############################################################################
// evolve_nse_plus_weak_rates().
//##############################################################################

void
evolve_nse_plus_weak_rates( nnt::Zone zone )
{

  gsl_vector * p_old_abundances, * p_abundance_changes;

  p_old_abundances = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  double d_ye_old = Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

  double d_log10_ye_old = log10( d_ye_old );

  double d_log10_ye =
    nnt::compute_1d_root(
      boost::bind(
        nse_plus_weak_rates_function,
        _1,
        boost::ref( zone ),
        d_ye_old
      ),
      d_log10_ye_old,
      1.1
    );

  nse_plus_weak_rates_function( d_log10_ye, zone, d_ye_old );
     
  p_abundance_changes =
    Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  gsl_vector_sub( p_abundance_changes, p_old_abundances );

  Libnucnet__Zone__updateAbundanceChanges(
    zone.getNucnetZone(),
    p_abundance_changes
  );

  gsl_vector_free( p_old_abundances );
  gsl_vector_free( p_abundance_changes );
  
}

//##############################################################################
// nse_plus_weak_rates_function().
//##############################################################################

double
nse_plus_weak_rates_function(
  double d_log10_ye,
  nnt::Zone& zone,
  double d_ye_old
)
{

  boost::tuple<double, double, double, double, double> T;

  double d_ye = pow( 10., d_log10_ye );

  zone.updateProperty(
    nnt::s_YE,
    d_ye
  ); 

  nnt::set_zone_abundances_to_equilibrium( zone );

  T = compute_all_yedot( d_ye, zone );

  return
    ( d_ye - d_ye_old ) /
      zone.getProperty<double>( nnt::s_DTIME )
    -
    T.get<4>();

}

//##############################################################################
// set_zone_for_evolution().
//##############################################################################

void
set_zone_for_evolution( nnt::Zone& zone )
{

  boost::any screening_data;
  boost::any coul_corr_data;

  if(
    Libnucnet__Zone__getScreeningFunction( zone.getNucnetZone() )
  )
  {
  
    screening_data =
      boost::any_cast<boost::function<boost::any()> >(
        zone.getFunction( nnt::s_SCREENING_DATA_FUNCTION )
      )();

    Libnucnet__Zone__setScreeningFunction(
      zone.getNucnetZone(),
      (Libnucnet__Zone__screeningFunction)
        Libnucnet__Zone__getScreeningFunction( zone.getNucnetZone() ),
      &screening_data
    );

  }
    
  if(
    Libnucnet__Zone__getNseCorrectionFactorFunction( zone.getNucnetZone() )
  )
  {
  
    coul_corr_data =
      boost::any_cast<boost::function<boost::any()> >(
        zone.getFunction( nnt::s_NSE_CORRECTION_FACTOR_DATA_FUNCTION )
      )();

    Libnucnet__Zone__setNseCorrectionFactorFunction(
      zone.getNucnetZone(),
      (Libnucnet__Species__nseCorrectionFactorFunction)
        Libnucnet__Zone__getNseCorrectionFactorFunction( zone.getNucnetZone() ),
      &coul_corr_data
    );

  }

  //--------------------------------------------------------------------------
  // Update data for reactions.
  //--------------------------------------------------------------------------

  if(
    zone.hasFunction( nnt::s_RATE_DATA_UPDATE_FUNCTION )
  )
  {

    boost::any_cast<boost::function<void( )> >(
      zone.getFunction( nnt::s_RATE_DATA_UPDATE_FUNCTION )
    )( );

  }

  //--------------------------------------------------------------------------
  // Compute rates.
  //--------------------------------------------------------------------------

  Libnucnet__Zone__computeRates(
    zone.getNucnetZone(),
    zone.getProperty<double>( nnt::s_T9 ),
    zone.getProperty<double>( nnt::s_RHO )
  ); 

  //--------------------------------------------------------------------------
  // Set weak detailed balance.
  //--------------------------------------------------------------------------

  set_weak_detailed_balance( zone );

  //--------------------------------------------------------------------------
  // Modify rates.
  //--------------------------------------------------------------------------

  if( zone.hasFunction( nnt::s_RATE_MODIFICATION_FUNCTION ) )
  {
    boost::any_cast<boost::function<void( )> >(
      zone.getFunction( nnt::s_RATE_MODIFICATION_FUNCTION )
    )( );
  }
  else
  {
    modify_rates( zone );
  }

  //--------------------------------------------------------------------------
  // Zero out small rates.
  //--------------------------------------------------------------------------

  if( zone.hasProperty( nnt::s_SMALL_RATES_THRESHOLD ) )
  {
    zero_out_small_rates(
      zone,
      zone.getProperty<double>( nnt::s_SMALL_RATES_THRESHOLD )
    );
  }

}

//##############################################################################
// get_evolution_matrix_and_vector().
//##############################################################################

std::pair< WnMatrix *, gsl_vector *>
get_evolution_matrix_and_vector( nnt::Zone& zone )
{

  set_zone_for_evolution( zone );

  //--------------------------------------------------------------------------
  // Return pair.
  //--------------------------------------------------------------------------

  return
    std::make_pair(
      Libnucnet__Zone__computeJacobianMatrix( zone.getNucnetZone() ),
      Libnucnet__Zone__computeFlowVector( zone.getNucnetZone() )
    );

}

//##############################################################################
// get_evolution_matrix().
//##############################################################################

WnMatrix *
get_evolution_matrix( nnt::Zone& zone )
{

  set_zone_for_evolution( zone );

  //--------------------------------------------------------------------------
  // Return matrix.
  //--------------------------------------------------------------------------

  return Libnucnet__Zone__computeJacobianMatrix( zone.getNucnetZone() );

}

//##############################################################################
// check_matrix_solution().
//##############################################################################

std::pair<double,double>
check_matrix_solution( nnt::Zone& zone, gsl_vector * p_sol )
{

  double d_checkT = 0, d_check = 0, d_abund, d_total = 0;
  size_t i_index;

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( zone.getNucnetZone() ) )
    );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_abund =
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      );

    i_index = Libnucnet__Species__getIndex( species.getNucnetSpecies() );

    if( zone.hasProperty( nnt::s_NEWTON_RAPHSON_ABUNDANCE ) )
    {
      if(
         d_abund >
         zone.getProperty<double>( nnt::s_NEWTON_RAPHSON_ABUNDANCE )
      )
      {
        d_checkT = fabs( gsl_vector_get( p_sol, i_index ) / d_abund );
        if( d_checkT > d_check ) d_check = d_checkT;
      }
    }
    else if( d_abund > D_Y_MIN )
    {
      d_checkT = fabs( gsl_vector_get( p_sol, i_index ) / d_abund );
      if( d_checkT > d_check ) d_check = d_checkT;
    }

    d_total +=
      gsl_pow_2(
        gsl_vector_get( p_sol, i_index ) *
        Libnucnet__Species__getA( species.getNucnetSpecies() )
      );
    
  }

  return std::make_pair( d_check, sqrt( d_total ) );

}
       
//##############################################################################
// network_t9_from_entropy_root(). 
//##############################################################################

double
network_t9_from_entropy_root(
  double d_x,
  nnt::Zone& zone
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

  safe_evolve(
    zone,
    zone.getProperty<double>( nnt::s_DTIME ),
    zone.getProperty<double>( nnt::s_DTIME )
  );

  d_result =
    compute_thermo_quantity(
      zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    )
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
    
//##############################################################################
// network_density_from_entropy_root(). 
//##############################################################################

double
network_density_from_entropy_root(
  double d_x,
  nnt::Zone& zone
)
{

  zone.updateProperty(
    nnt::s_RHO,
    d_x
  );

  return
    compute_thermo_quantity(
      zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    )
    -      
    zone.getProperty<double>( nnt::s_ENTROPY_PER_NUCLEON );

}

} // namespace user
