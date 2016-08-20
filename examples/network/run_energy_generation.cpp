////////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code for running a network calculation with energy
//!        generation.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <algorithm>

#include <Libnucnet.h>

#include <boost/bind.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/program_options.hpp>

#include <boost/format.hpp>
#include <boost/assign.hpp>

#include "nnt/two_d_weak_rates.h"
#include "user/remove_duplicate.h"
#include "user/user_rate_functions.h"
#include "user/network_limiter.h"
#include "user/flow_utilities.h"

#include "user/hydro_helper.h"

typedef user::state_type my_state_type;

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_DT0          1.e-15  /* Initial time step */
#define D_REG_T        0.15    /* Time step change regulator for dt update */
#define D_REG_Y        0.15    /* Abundance change regulator for dt update */
#define D_X_REG_T      0.15    /* x change regulator for dt update */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */
#define D_LIM_CUTOFF   1.e-50  /* Cutoff abundance for network limiter */

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

#define VALIDATE       "no"

//##############################################################################
// Strings.
//##############################################################################

#define S_SOLVER       nnt::s_ARROW // Solver type: ARROW or GSL

#define S_P_0          "P_0"

#define S_X            "x"

#define S_PARTICLE     nnt::s_TOTAL

#define S_OBSERVE  "observe"

#define S_ACCELERATION_FUNCTION  "acceleration function"
#define S_ENERGY_GENERATION_FUNCTION  "energy generation function"
#define S_EVOLVE_FUNCTION "evolution function"
#define S_OBSERVER_FUNCTION  "observer function"
#define S_RHO_FUNCTION  "rho function"
#define S_D_LN_RHO_DT_FUNCTION "d ln rho dt function"

namespace po = boost::program_options;

//##############################################################################
// energy_generation_rhs.
//##############################################################################

class energy_generation_rhs
{

  nnt::Zone& zone;
  Libnucnet__NetView * pView;

  public:
    energy_generation_rhs( nnt::Zone& _zone, Libnucnet__NetView * p_view ) :
      zone( _zone ), pView( p_view ) {}

    void
    operator() ( const my_state_type &x, my_state_type &dxdt, const double d_t )
    {

      double d_dt, d_energy_generation; //, d_energy_loss;
      gsl_vector * p_abundances, * p_abundance_changes;

      p_abundances = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

      p_abundance_changes =
        Libnucnet__Zone__getAbundanceChanges( zone.getNucnetZone() );

      d_dt =
        d_t 
        -
        zone.getProperty<double>( nnt::s_TIME )
        + 1.e-25;

      zone.updateProperty(
        nnt::s_DTIME,
        d_dt
      );

      zone.updateProperty(
        nnt::s_T9,
        x[2]
      );

      zone.updateProperty(
        nnt::s_RHO,
        boost::any_cast<boost::function<double( const my_state_type& )> >(
          zone.getFunction( S_RHO_FUNCTION )
        )( x )
      );

      double d_ln_rho_dt =
        boost::any_cast<boost::function<double( const my_state_type& )> >(
          zone.getFunction( S_D_LN_RHO_DT_FUNCTION )
        )( x );

      boost::any_cast<
        boost::function<void( Libnucnet__NetView *, const double )>
      >(
        zone.getFunction( S_EVOLVE_FUNCTION )
      )( pView, d_dt );

      dxdt[0] = x[1];

      dxdt[1] = 
        boost::any_cast<
          boost::function<double( const my_state_type&, const double )>
        >(
          zone.getFunction( S_ACCELERATION_FUNCTION )
        )( x, d_t );

      d_energy_generation =
        boost::any_cast<boost::function<double( Libnucnet__NetView * )> >(
          zone.getFunction( S_ENERGY_GENERATION_FUNCTION )
        )( pView );

/*
      d_energy_loss =
        boost::any_cast<boost::function<double(
          const state_type&,
          const double
        )> >(
          zone.getFunction( S_ENERGY_LOSS_FUNCTION )
        )( x, d_t );
*/

        dxdt[2] =
          ( 
            1. / 
            (
              user::compute_thermo_quantity(
                zone,
                nnt::s_SPECIFIC_HEAT_PER_NUCLEON, 
                nnt::s_TOTAL
              ) * GSL_CONST_CGSM_BOLTZMANN
            ) 
          )
          *
          (
            zone.getProperty<double>( nnt::s_T9 ) *
               GSL_CONST_NUM_GIGA
               * user::compute_thermo_quantity(
                   zone, nnt::s_DPDT, nnt::s_TOTAL
                 )
               * d_ln_rho_dt
            / 
            (
              zone.getProperty<double>( nnt::s_RHO ) * GSL_CONST_NUM_AVOGADRO
            )
            + d_energy_generation // - d_energy_loss
          );

        dxdt[2] /= GSL_CONST_NUM_GIGA;

        if( zone.hasFunction( S_OBSERVER_FUNCTION ) )
        {
          boost::any_cast<
            boost::function<
              void(
                const my_state_type&,
                const my_state_type&,
                const double
              )
            >
          >( zone.getFunction( S_OBSERVER_FUNCTION ) )( x, dxdt, d_t );
        }

        Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abundances );
      
        Libnucnet__Zone__updateAbundanceChanges(
          zone.getNucnetZone(),
          p_abundance_changes
        );

    }
      
}; 

//##############################################################################
// get_input().
//##############################################################################

std::pair<Libnucnet *, Libnucnet__NetView *>
get_input( int argc, char **argv )
{

  Libnucnet * p_nucnet;
  Libnucnet__NetView * p_view = NULL;

  try
  {

    std::string s_nuc_xpath = "", s_reac_xpath = "";
    std::string s_enuc_nuc_xpath = "", s_enuc_reac_xpath = "";

    std::string s_purpose = "\nPurpose: run a network calculation with entropy generation for the input xml_file for the selected nuclei and reactions and for the selected nuclei and reactions for entropy generation.";

    po::options_description desc("\nAllowed options");
    desc.add_options()
      ( "help", "print out this help message and exit" )
      (
       "nuc_xpath",
       po::value<std::string>(),
       "XPath to select nuclides (default: all nuclides)"
      )
      (
       "reac_xpath",
       po::value<std::string>(),
       "XPath to select reaction (default: all reactions)"
      )
      (
       "enuc_nuc_xpath",
       po::value<std::string>(),
       "XPath to select nuclides for entropy generation (default: a step's evolution network nuclides)"
      )
      (
       "enuc_reac_xpath",
       po::value<std::string>(),
       "XPath to select reactions for entropy generation (default: a step's evolution network reactions)"
      )
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line( argc, argv, desc), vm );
    po::notify(vm);    

    if( argc < 4 || vm.count("help") == 1 )
    {
      std::cerr <<
        "\nUsage: " << argv[0] << " net_xml zone_xml output_xml [options]" <<
        std::endl;
      std::cerr << s_purpose << std::endl;
      std::cout << desc << "\n";
      exit( EXIT_FAILURE );
    }

    if( vm.count("nuc_xpath") == 1 )
    {
      s_nuc_xpath = vm["nuc_xpath"].as<std::string>();
    }

    if( vm.count("reac_xpath") == 1 )
    {
      s_reac_xpath = vm["reac_xpath"].as<std::string>();
    }

    if( vm.count("enuc_nuc_xpath") == 1 )
    {
      s_enuc_nuc_xpath = vm["enuc_nuc_xpath"].as<std::string>();
    }

    if( vm.count("enuc_reac_xpath") == 1 )
    {
      s_enuc_reac_xpath = vm["enuc_reac_xpath"].as<std::string>();
    }

    //==========================================================================
    // Validate input file.
    //==========================================================================

    if( strcmp( VALIDATE, "yes" ) == 0 )
    {
      if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
        fprintf( stderr, "Not valid libnucnet input!\n" );
        exit( EXIT_FAILURE );
      }
    }

    //==========================================================================
    // Get network and view.
    //==========================================================================

    p_nucnet = Libnucnet__new();

    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      s_nuc_xpath.c_str(),
      s_reac_xpath.c_str()
    );

    Libnucnet__assignZoneDataFromXml(
      p_nucnet,
      argv[2],
      ""
    );

    if( vm.count("enuc_nuc_xpath") == 1 || vm.count("enuc_reac_xpath") == 1 )
    {
      p_view =
        Libnucnet__NetView__new(
          Libnucnet__getNet( p_nucnet ),
          s_enuc_nuc_xpath.c_str(),
          s_enuc_reac_xpath.c_str()
        );
    }

    return std::make_pair( p_nucnet, p_view );

  }
  catch( std::exception& e )
  {
    std::cerr << "error: " << e.what() << "\n";
    exit( EXIT_FAILURE );
  }
  catch(...)
  {
    std::cerr << "Exception of unknown type!\n";
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step, k = 0;
  double d_t, d_dt;
  Libnucnet * p_my_nucnet, * p_my_output;
  Libnucnet__NetView * p_view;
  nnt::Zone zone;
  char s_property[32];
  std::set<std::string> isolated_species_set;

  my_state_type
    x(3), xold(3), x_lim = boost::assign::list_of(1.e-10)(1.)(1.e-5);

  //============================================================================
  // Check input.
  //============================================================================

  boost::tie( p_my_nucnet, p_view ) = get_input( argc, argv );

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Set the zone.
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  //============================================================================
  // Use approximate weak rates or not.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) &&
    zone.getProperty<std::string>( nnt::s_USE_APPROXIMATE_WEAK_RATES ) == "yes"   )
  {
    user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );
  }

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

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
  // Remove isolated species if desired.
  //============================================================================

  isolated_species_set =
    user::get_isolated_species(
      Libnucnet__getNet( p_my_nucnet ),
      "",
      ""
    );

  BOOST_FOREACH( std::string s_species, isolated_species_set )
  {

//  Careful that you don't remove a species with non-zero abundance!

    std::cout << s_species << std::endl;

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
  // Set the acceleration function.
  // Must return the scaled acceleration for the given available data.
  //============================================================================

  zone.updateFunction(
    S_ACCELERATION_FUNCTION,
    static_cast<boost::function<double( const my_state_type&, const double )> >(
      boost::bind(
        user::acceleration,
        boost::ref( zone ),
        _1,
        _2
      )
    )
  );

  //============================================================================
  // Set the rho function.
  // Must return the density for the given available data.
  //============================================================================

  zone.updateFunction(
    S_RHO_FUNCTION,
    static_cast<boost::function<double( const my_state_type& )> >(
      boost::bind(
        user::rho_function,
        boost::ref( zone ),
        _1
      )
    )
  );

  //============================================================================
  // Set the d(ln rho)/dt function.
  // Must return the time rate of change of the ln of the density for the
  // given available data.
  //============================================================================

  zone.updateFunction(
    S_D_LN_RHO_DT_FUNCTION,
    static_cast<boost::function<double( const my_state_type& )> >(
      boost::bind(
        user::d_ln_rho_dt_function,
        boost::ref( zone ),
        _1
      )
    )
  );

  //============================================================================
  // Set the energy generation function.
  //============================================================================

  zone.updateFunction(
    S_ENERGY_GENERATION_FUNCTION,
    static_cast<boost::function<double( Libnucnet__NetView * )> >(
      boost::bind(
        user::compute_energy_generation_rate_per_nucleon,
        boost::ref( zone ),
        _1
      )
    )
  );

  //============================================================================
  // Set the abundance evolver.
  // Must evolve the abundances over a given input timestep.
  //============================================================================

  zone.updateFunction(
    S_EVOLVE_FUNCTION,
    static_cast<boost::function<void( Libnucnet__NetView *, const double )> >(
      boost::bind(
        user::evolve_function,
        boost::ref( zone ),
        _1,
        _2
      )
    )
  );

  //============================================================================
  // Set the observer function with basic prototype
  //   void( const my_state_type& x, const my_state_type& dxdt, const double t )
  //============================================================================

  if( zone.hasProperty( S_OBSERVE ) && 
      zone.getProperty<std::string>( S_OBSERVE ) == "yes"
  )
  {
    zone.updateFunction(
      S_OBSERVER_FUNCTION,
      static_cast<
        boost::function<
          void(
            const my_state_type&,
            const my_state_type&,
            const double
          )
        >
      >( boost::bind(
           user::observer_function,
           boost::ref( zone ),
           _1,
           _2,
           _3
         )
      )
    );
  }
            
  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( S_SOLVER == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    zone.updateProperty( nnt::s_SOLVER, nnt::s_ARROW );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "3" );

  }

  //============================================================================
  // Create output.
  //============================================================================

  p_my_output = nnt::create_network_copy( p_my_nucnet );

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_output,
    "%.15e"
  );

  //============================================================================
  // Set initial temperature and density.
  //============================================================================

  zone.updateProperty(
    nnt::s_T9,
    zone.getProperty<std::string>( nnt::s_T9_0 )
  );
  
  zone.updateProperty(
    nnt::s_RHO,
    zone.getProperty<std::string>( nnt::s_RHO_0 )
  );
  
  //============================================================================
  // Initialize the system.
  //============================================================================

  zone.updateProperty( nnt::s_PARTICLE, S_PARTICLE);

  if( zone.hasProperty( nnt::s_DTIME ) )
    d_dt = zone.getProperty<double>( nnt::s_DTIME );
  else
    d_dt = D_DT0;

  d_t = 0.0;

  x[0] = zone.getProperty<double>( S_X, "0" );
  x[1] = x[0] / ( 3. * zone.getProperty<double>( nnt::s_TAU ) );
  x[2] = zone.getProperty<double>( nnt::s_T9 );

  zone.updateProperty(
    S_P_0,
    user::compute_thermo_quantity(
      zone,
      nnt::s_PRESSURE,
      zone.getProperty<std::string>( nnt::s_PARTICLE )
    )
  );

  user::limit_evolution_network( zone, D_LIM_CUTOFF );

  //============================================================================
  // Choose the stepper.
  //============================================================================

  boost::numeric::odeint::adams_bashforth< 8, my_state_type > stepper;

  //============================================================================
  // Evolve network while t < final t. 
  //============================================================================

  i_step = 0;

  while ( d_t < zone.getProperty<double>( nnt::s_TEND ) )
  {

  //============================================================================
  // Set time.
  //============================================================================

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );  

  //============================================================================
  // Save old values.
  //============================================================================

    std::copy( x.begin(), x.end(), xold.begin() );

  //============================================================================
  // Evolve step.
  //============================================================================

    Libnucnet__NetView * p_enuc_view;

    if( p_view )
      p_enuc_view = p_view;
    else
      p_enuc_view = zone.getNetView( EVOLUTION_NETWORK );

    energy_generation_rhs my_rhs( zone, p_enuc_view );

    stepper.do_step( my_rhs, x, d_t, d_dt );

  //============================================================================
  // Update properties.
  //============================================================================

    d_t += d_dt;

    zone.updateProperty( nnt::s_DTIME, d_dt );

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

    zone.updateProperty(
      nnt::s_RHO,
      user::rho_function( zone, x )
    );

    zone.updateProperty(
      nnt::s_T9,
      x[2]
    );

    boost::any_cast<
      boost::function<void( Libnucnet__NetView *, const double )>
    >(
      zone.getFunction( S_EVOLVE_FUNCTION )
    )( zone.getNetView( EVOLUTION_NETWORK ), d_dt );

  //============================================================================
  // Output step data.
  //============================================================================

    if(
      zone.hasProperty( S_OBSERVE ) &&
      zone.getProperty<std::string>( S_OBSERVE ) == "yes"
    )
    {
      std::cout <<
        boost::format( "t = %g, x = {%g, %g, %g}\n\n" ) %
        d_t % x[0] % x[1] % x[2];
      std::cout << boost::format( "-----------\n\n" );
    }

    zone.updateProperty( S_X, "0", x[0] );

    zone.updateProperty( S_X, "1", x[1] );

  //============================================================================
  // Print out abundances.
  //============================================================================

    if( ( i_step++ % zone.getProperty<int>( nnt::s_STEPS ) ) == 0 ||
        d_t >= zone.getProperty<double>( nnt::s_TEND )
    )
    {
      sprintf( s_property, "%d", ++k );
      Libnucnet__relabelZone(
        p_my_nucnet,
        zone.getNucnetZone(),
        s_property,
        NULL,
        NULL
      );
      nnt::print_zone_abundances( zone );
      nnt::write_xml( p_my_output, zone.getNucnetZone() );
      Libnucnet__writeToXmlFile( p_my_output, argv[3] );
    }

  //============================================================================
  // Limit network.
  //============================================================================

  user::limit_evolution_network( zone, D_LIM_CUTOFF );

  //============================================================================
  // Update timestep.
  //============================================================================

    double d_h = 1.e99;
    for( size_t i = 0; i < x.size(); i++ )
    {
      double delta = fabs( ( x[i] - xold[i] ) / x[i] );
      if( delta > 0 && fabs( x[i] ) > x_lim[i] )
        d_h = GSL_MIN( d_h, D_X_REG_T * d_dt / delta );
    }

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if( d_dt > d_h ) d_dt = d_h;

    if ( d_t + d_dt > zone.getProperty<double>( nnt::s_TEND ) ) {

      d_dt = zone.getProperty<double>( nnt::s_TEND ) - d_t;

    }

  }  

  //============================================================================
  // Write output.
  //============================================================================

  Libnucnet__writeToXmlFile( p_my_output, argv[3] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
