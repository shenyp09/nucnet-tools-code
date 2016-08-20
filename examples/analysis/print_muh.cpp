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
//! \brief Example code to print out the heavy nuclide chemical potential
//!     in zones from a network xml file.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include <Libnuceq.h>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

#include "user/network_utilities.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output.xml \"[position() >= last() - 10]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " computes quantities related to heavy nuclei for selected zones" <<
      std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones of interest\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );

  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnuceq * p_my_equil;
  Libnuceq__Cluster * p_my_cluster;
  const char S_XPATH[] = "[z > 2]";

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
      argv[2]
    );

  //============================================================================
  // Create equilibrium and cluster.
  //============================================================================

  p_my_equil =
    Libnuceq__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  p_my_cluster =
    Libnuceq__newCluster(
      p_my_equil,
      S_XPATH
    );

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t my_zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, my_zone_list )
  {

    //--------------------------------------------------------------------------
    // Set the Ye for the equilibrium.
    //--------------------------------------------------------------------------

    Libnuceq__setYe(
      p_my_equil,
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
    );

    //--------------------------------------------------------------------------
    // Compute yh and update the cluster constraint.
    //--------------------------------------------------------------------------

    Libnuceq__Cluster__updateConstraint(
      p_my_cluster,
      user::compute_cluster_abundance_moment( zone, S_XPATH, "a", 0 )
    );

    //--------------------------------------------------------------------------
    // Compute equilibrium and print out.
    //--------------------------------------------------------------------------

    Libnuceq__computeEquilibrium(
      p_my_equil,
      zone.getProperty<double>( nnt::s_T9 ),
      zone.getProperty<double>( nnt::s_RHO )
    );

    std::cout <<
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) << "  " <<
      zone.getProperty<std::string>( nnt::s_TIME ) << "  " <<
      zone.getProperty<std::string>( nnt::s_T9 ) << "  " <<
      zone.getProperty<std::string>( nnt::s_RHO ) << "  " <<
      Libnuceq__Cluster__getConstraint( p_my_cluster ) << "  " <<
      Libnuceq__Cluster__getMukT( p_my_cluster ) <<
      std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnuceq__free( p_my_equil );

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
