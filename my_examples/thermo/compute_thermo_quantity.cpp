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
//! \brief Example code to compute thermodynamic quantities from a network xml
//!        file.
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <sstream>
#include <vector>
#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/thermo.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--help" ) == 0 )
  {
    nnt::Zone zone;
    std::cerr << "The quantities available to compute are:" << std::endl;
    std::cerr << std::endl;
    std::vector<std::string> quantities =
      user::list_zone_thermo_quantities( zone );
    BOOST_FOREACH( std::string s, quantities )
    {
      std::cerr << "  " << s << std::endl;
    }
    std::cerr << std::endl;
    std::cerr << "For more information about a specific quantity, type, for example," <<
      std::endl << std::endl;
    std::cerr << argv[0] << " --help \"baryon pressure\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc == 3 && strcmp( argv[1], "--help" ) == 0 )
  {
    nnt::Zone zone;
    std::cerr << std::endl;
    {
      std::cerr <<
        user::get_zone_thermo_quantity_doc( zone, argv[2] ) << std::endl;
    }
    std::cerr << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cerr << std::endl;
    std::cerr << argv[0] << " my_output.xml \"baryon entropy per nucleon\" \"electron internal energy density\" \"photon pressure\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc <= 2 ) {
    std::cerr << std::endl;
    std::cerr << "Purpose: " << argv[0] <<
      " computes thermodynamic quantities from a network output xml file" <<
      std::endl;

    fprintf(
      stderr, "\nUsage: %s file quantity ...\n\n", argv[0]
    );
    fprintf(
      stderr, "  file = input xml filename\n\n"
    );
    fprintf(
      stderr, "  quantity = quantity to compute (as many as desired)\n\n"
    );
    std::cerr << "For an example usage, type " << std::endl << std::endl;
    std::cerr << argv[0] << " --example" << std::endl << std::endl;
    std::cerr << "For the available quantities to compute, type " <<
      std::endl << std::endl;
    std::cerr << argv[0] << " --help" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char * argv[] ) {

  std::list<std::string> quantity_list;
  Libnucnet * p_my_nucnet;

  //============================================================================
  // Check input.
  //============================================================================
  
  check_input( argc, argv ); 

  //============================================================================
  // Read input data.
  //============================================================================

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      NULL
    );

  //============================================================================
  // Set quantities.
  //============================================================================

  for( int i = 2; i < argc; i++ )
    quantity_list.push_back( argv[i] );

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
    std::cout << Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 );
    BOOST_FOREACH( std::string s_quantity, quantity_list )
    {
      std::cout << "  " <<
        user::compute_thermo_quantity( zone, s_quantity );
    }
    std::cout << std::endl;
  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
