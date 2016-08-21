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
//! \brief Example code to convergence of solution for one time step.
////////////////////////////////////////////////////////////////////////////////

/*##############################################################################
// Includes.
//############################################################################*/

#include <boost/program_options.hpp>

#include <Libnucnet.h>

#include "user/evolve.h"
#include "user/remove_duplicate.h"
#include "user/matrix_solver.h"
#include "user/user_rate_functions.h"

namespace po = boost::program_options;

/*##############################################################################
// Enumeration.
//############################################################################*/

enum solvers { ARROW, GSL };

/*##############################################################################
// Screening.  "no" = no screening, "yes" = screening.
//############################################################################*/

#define VALIDATE       "no"

/*##############################################################################
// Prototypes.
//############################################################################*/

int sort_function( Libnucnet__Species *, Libnucnet__Species * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  nnt::Zone zone;
  WnMatrix *p_matrix;
  gsl_vector *p_rhs, *p_sol, *p_y_old, *p_work;
  int i_continue = 1;
  std::string s_nuc_xpath = "", s_reac_xpath = "";


  //============================================================================
  // Check input.
  //============================================================================

  try
  {

    std::string s_purpose = "\nPurpose: illustrate the convergence of the matrix solution for the input net_xml and zone_xml files for the selected nuclei and reactions.";

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
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line( argc, argv, desc), vm );
    po::notify(vm);    

    if( argc < 3 || vm.count("help") == 1 )
    {
      std::cout <<
        "\nUsage: " << argv[0] << " net_xml zone_xml [options]" << std::endl;
      std::cout << s_purpose << std::endl;
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

  //============================================================================
  // Validate input file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_input_xml( argv[1] ) )
    {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      return EXIT_FAILURE;
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Net__updateFromXml( 
    Libnucnet__getNet( p_my_nucnet ),
    argv[1],
    s_nuc_xpath.c_str(),
    s_reac_xpath.c_str()
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  //============================================================================
  // Register user-supplied rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions(
    Libnucnet__getNet( p_my_nucnet )
  );

  //============================================================================
  // Get the zone.
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  //============================================================================
  // Save the abundances.
  //============================================================================

  p_y_old = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  //============================================================================
  // Newton-Raphson Iterations.
  //============================================================================

  while( i_continue )
  {

    //--------------------------------------------------------------------------
    // Get matrix and rhs vector.
    //--------------------------------------------------------------------------

    boost::tie( p_matrix, p_rhs ) =
      user::get_evolution_matrix_and_vector( zone );

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
    gsl_vector_scale( p_work, 1. / zone.getProperty<double>( nnt::s_DTIME ) );
    gsl_vector_sub( p_rhs, p_work );

    gsl_vector_free( p_work );

    //--------------------------------------------------------------------------
    // Solve matrix equation.
    //--------------------------------------------------------------------------

    p_sol = user::solve_matrix_for_zone( zone, p_matrix, p_rhs );

    //==========================================================================
    // Solve and print out.
    //==========================================================================

    p_work = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

    fprintf(
      stdout,
      "\n    i   Z   A   Y_old[i]          dY[i]        Y_new[i]\n"
    );

    BOOST_FOREACH(
      nnt::Species species,
      nnt::make_species_list(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    )
    {

      size_t i = Libnucnet__Species__getIndex( species.getNucnetSpecies() );

      fprintf(
        stdout,
        "%5lu %5u %5u   %.5e    %.5e    %.5e\n",
        (unsigned long) i,
        Libnucnet__Species__getZ( species.getNucnetSpecies() ),
        Libnucnet__Species__getA( species.getNucnetSpecies() ),
        gsl_vector_get( p_work, i ),
        gsl_vector_get( p_sol, i ),
        gsl_vector_get( p_work, i ) + gsl_vector_get( p_sol, i )
      );

    }

    //--------------------------------------------------------------------------
    // Update abundances.
    //--------------------------------------------------------------------------

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

    std::cout << std::endl << "quit [0]   continue [1]" << std::endl;
    std::cin >> i_continue;

  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_y_old );

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
