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
//! \brief Example code for running a network calculation at constant entropy.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>
#include "nnt/auxiliary.h"
#include "nnt/math.h"
#include "nnt/string_defs.h"
#include "user/evolve.h"
#include "user/hydro.h"
#include "user/thermo.h"
#include "user/user_rate_functions.h"
#include "user/neutrino_rate_functions.h"
#include "user/network_utilities.h"
#include "user/remove_duplicate.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define MY_BUF_SIZE    32      /* String buffer size */
#define D_DT0          1.e-05  /* Initial time step */
#define D_REG_T        0.15    /* Time step change regulator for dt update */
#define D_REG_Y        0.15    /* Abundance change regulator for dt update */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */

//##############################################################################
// Strings.
//##############################################################################

#define S_SOLVER_TYPE   nnt::s_ARROW  /* nnt::s_ARROW or nnt::s_GSL */
#define S_T9_PREVIOUS   "t9 prev"
#define S_SPECIES_REMOVAL_NUC_XPATH  "species removal nuclide xpath"
#define S_SPECIES_REMOVAL_REAC_XPATH  "species removal reaction xpath"

//##############################################################################
// Booleans.
//##############################################################################

#define B_OUTPUT_EVERY_TIME_DUMP    false  // Change to true to write to xml
                                           // every time dump.  False just
                                           // writes output at end of
                                           // calculation.

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step;
  double d_tend, d_t, d_dt;
  Libnucnet *p_my_nucnet, *p_my_output;

  nnt::Zone zone;

  //============================================================================
  // Get the nucnet.
  //============================================================================
  
  p_my_nucnet = user::get_nucnet( argc, argv );
  
  //============================================================================
  // Get the zone.
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  if( !zone.getNucnetZone() )
  {
    fprintf( stderr, "No input zone with labels (0,0,0).\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Register user-supplied rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  user::register_neutrino_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

  //============================================================================
  // Update with approximate weak rates if desired.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) &&
    zone.getProperty<std::string>( nnt::s_USE_APPROXIMATE_WEAK_RATES ) == "yes"   )
  {
      user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );
  }

  //============================================================================
  // If not set, set the neutrino chemical potential to -inf.
  //============================================================================
  
  if( !zone.hasProperty( nnt::s_MU_NUE_KT ) )
  {
    zone.updateProperty(
      nnt::s_MU_NUE_KT,
      "-inf"
    );
  }

  //============================================================================
  // Set solver.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(), nnt::s_SOLVER, NULL, NULL, S_SOLVER_TYPE
  );

  //============================================================================
  // Remove isolated species since we start at NSE.  NSE may set isolated
  // species to have non-zero abundance which will then not evolve.
  //============================================================================

  std::set<std::string> isolated_species_set =
    user::get_isolated_species(
      Libnucnet__getNet( p_my_nucnet ),
      "",
      ""
    );

  BOOST_FOREACH( std::string s_species, isolated_species_set )
  {

    std::cout << "Removing isolated species: " << s_species << std::endl;

    Libnucnet__Nuc__removeSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc(
          Libnucnet__getNet( p_my_nucnet )
        ),
        s_species.c_str()
      )
    );

  }
    
  //============================================================================
  // Set screening, Coulomb correction, and rate update functions.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_USE_SCREENING ) &&
    zone.getProperty<std::string>( nnt::s_USE_SCREENING ) == "yes" 
  )
  {
    user::set_screening_function( zone );
  }

  if(
    zone.hasProperty( nnt::s_USE_NSE_CORRECTION ) &&
    zone.getProperty<std::string>( nnt::s_USE_NSE_CORRECTION ) == "yes" 
  )
  {
    user::set_nse_correction_function( zone );
  }

  user::set_rate_data_update_function( zone );

  //============================================================================
  // Create output network structure.
  //============================================================================

  p_my_output = nnt::create_network_copy( p_my_nucnet );

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( zone.getProperty<std::string>( nnt::s_SOLVER ) == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function)
         nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_output ) ),
      (Libnucnet__Species__compare_function)
         nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_output ) )
    );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "3" );

  }

  //============================================================================
  // Initialize temperature parameters.
  //============================================================================

  zone.updateProperty(
    nnt::s_T9_0,
    zone.getProperty<std::string>( nnt::s_T9 )
  );

  zone.updateProperty(
    S_T9_PREVIOUS,
    zone.getProperty<std::string>( nnt::s_T9 )
  );

  //============================================================================
  // Start with network in NSE at input t9 and Ye.
  //============================================================================

  zone.updateProperty(
    nnt::s_RHO,
    pow(
      10.,
      nnt::compute_1d_root(
        boost::bind(
          user::compute_log10_density_entropy_root_with_equilibrium,
          _1,
          boost::ref( zone )
        ),
        log10( zone.getProperty<double>( nnt::s_T9 ) ),
        1.1
      )
    )
  );

  zone.updateProperty(
    nnt::s_RHO_0,
    zone.getProperty<std::string>( nnt::s_RHO )
  );

  //============================================================================
  // Limit network.
  //============================================================================

  user::limit_evolution_network( zone );

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;
  d_t = zone.getProperty<double>( nnt::s_TIME );

  d_tend = zone.getProperty<double>( nnt::s_TEND );

  while( d_t < d_tend )
  {

    d_dt = zone.getProperty<double>( nnt::s_DTIME );

    d_t += d_dt;

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

  //============================================================================
  // Update density and radius.  Swap neutrinos if below resonance.
  //============================================================================

    user::update_zone_properties( zone );

    if( zone.hasProperty( S_RHO_RES ) )
      user::swap_neutrinos( zone );

  //============================================================================
  // Evolve abundances.
  //============================================================================

    zone.updateProperty(
      nnt::s_T9,
      nnt::compute_1d_root(
        boost::bind(
          user::network_t9_from_entropy_root,
          _1,
          boost::ref( zone )
        ),
        zone.getProperty<double>( S_T9_PREVIOUS ),
        1.001
      )
    );

    user::evolve( zone );

    zone.updateProperty(
      S_T9_PREVIOUS,
      zone.getProperty<std::string>( nnt::s_T9 )
    );

  //============================================================================
  // Print out abundances.
  //============================================================================

    if(
      (
        i_step % zone.getProperty<int>( nnt::s_STEPS) == 0 
      ) ||
        d_t >= d_tend
    )
    {
      zone.updateProperty(
        nnt::s_YE,
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
      );
      nnt::print_zone_abundances(zone );
      nnt::write_xml( p_my_output, zone.getNucnetZone() );
      if( B_OUTPUT_EVERY_TIME_DUMP )
      {
        Libnucnet__writeToXmlFile( p_my_output, argv[3] );
      }
    }

  //============================================================================
  // Update timestep.
  //============================================================================

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if ( d_t + d_dt > d_tend ) {

      d_dt = d_tend - d_t;

    }

    zone.updateProperty(
      nnt::s_DTIME,
      d_dt
    );

    user::limit_evolution_network( zone );

    nnt::normalize_zone_abundances( zone );

    i_step++;

  }  

  //============================================================================
  // Write out and free p_my_output.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );
  Libnucnet__writeToXmlFile( p_my_output, argv[3] );
  Libnucnet__free( p_my_output );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
