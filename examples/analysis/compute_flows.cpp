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
//! \brief Example code to compute the flows for given zones (chosen by
//!    XPath) in a network xml file.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"
#include "nnt/two_d_weak_rates.h"

#include "user/network_utilities.h"
#include "user/flow_utilities.h"
#include "user/weak_utilities.h"
#include "user/user_rate_functions.h"
#include "user/neutrino_rate_functions.h"

#ifdef MY_USER
#include "my_user/my_rates.h"
#endif

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__NetView * p_net_view;
  std::pair<double,double> flows;
  double d_min_flow;

  //============================================================================
  // Check input.
  //============================================================================

   if ( argc < 3 || argc > 5 ) {
      fprintf(
        stderr, "\nUsage: %s xml_filename zone_xpath reac_xpath min_flow\n\n", argv[0]
      );
      fprintf(
        stderr, "  xml_filename = input xml filename\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = XPath to select zones for flows\n\n"
      );
      fprintf(
        stderr, "  reac_xpath = reaction xpath (optional--required if min_flow present)\n\n"
      );
      fprintf(
        stderr, "  min_flow = minimum net flow to print out (optional--if not set, min is 1.e-50)\n\n"
      );
      return EXIT_FAILURE;
   }

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[2] );

  if( !p_my_nucnet ) {
    fprintf( stderr, "Input data not read!\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac(
      Libnucnet__getNet( p_my_nucnet )
    )
  );

  user::register_neutrino_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

#ifdef MY_USER
  my_user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );
#endif

  //============================================================================
  // Get the zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  //============================================================================
  // Get valid net view.
  //============================================================================

  if( argc == 3 )
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", "" );
  else
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", argv[3] );

  //============================================================================
  // Set minimum flow.
  //============================================================================

  if( argc == 5 )
    d_min_flow = boost::lexical_cast<double>( argv[4] );
  else
    d_min_flow = 1.e-50;

  //============================================================================
  // Set compare function.
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) ),
    (Libnucnet__Reaction__compare_function)
       nnt::compare_reactions_by_string
  );

  //============================================================================
  // Get the reactions.
  //============================================================================

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );

  Libnucnet__NetView__free( p_net_view );

  //============================================================================
  // Iterate the zones.
  //============================================================================

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
   
    if( !( zone == *(zone_list.begin()) ) )
    {
      Libnucnet__Zone__copy_net_views(
        zone.getNucnetZone(),
        (*zone_list.begin()).getNucnetZone()
      );
    }

    user::update_rate_functions_data( zone );

    if(
      zone.hasProperty( nnt::s_USE_SCREENING )
      &&
      zone.getProperty<std::string>( nnt::s_USE_SCREENING ) == "yes"
    )
    {
      user::set_screening_function( zone );
    }

    if(
      zone.hasProperty( nnt::s_USE_NSE_CORRECTION )
      &&
      zone.getProperty<std::string>( nnt::s_USE_NSE_CORRECTION ) == "yes"
    )
    {
      user::set_nse_correction_function( zone );
    }

    user::flow_data_tuple_t my_data_tuple = user::make_flow_data_tuple( zone );

    //==========================================================================
    // Print conditions.
    //==========================================================================

    if( zone.hasProperty( nnt::s_TIME ) )
      std::cout <<
        "time(s) = " << zone.getProperty<std::string>( nnt::s_TIME ) << " " <<
        "t9 = " << zone.getProperty<std::string>( nnt::s_T9 ) << " " <<
        "rho(g/cc) = " << zone.getProperty<std::string>( nnt::s_RHO ) << " " <<
        std::endl;
    else
      std::cout <<
        "t9 = " << zone.getProperty<std::string>( nnt::s_T9 ) << " " <<
        "rho(g/cc) = " << zone.getProperty<std::string>( nnt::s_RHO ) << " " <<
        std::endl;

    fprintf( stdout, "\n\t\t\tReaction\t\t\t  Forward     Reverse     Net\n" );
    printf(
      "=======================================================   =========   =========   =========\n"
    );

    //==========================================================================
    // Print flows.
    //==========================================================================
  
    BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
    {

      flows =
        user::compute_flows_for_reaction(
          zone,
          reaction.getNucnetReaction(),
          my_data_tuple
      );

      if( fabs( flows.first - flows.second ) >= d_min_flow )
      {

        fprintf(
      	  stdout,
	  "%-55s%12.3e%12.3e%12.3e\n",
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
          flows.first,
          flows.second,
          flows.first - flows.second
        );

      }

    }

    fprintf( stdout, "\n" );

  }

  //===========================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
