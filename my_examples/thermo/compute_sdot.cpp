//////////////////////////////////////////////////////////////////////////////
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
//!
//! \file
//! \brief Example code to compute the rate of change of the entropy per
//!    nucleon in zones in a output network xml file.
//!
////////////////////////////////////////////////////////////////////////////////

#include <boost/program_options.hpp>

#include <Libnucnet.h>

#include "nnt/wrappers.hpp"
#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

#include "user/aa522a25.h"
#include "user/flow_utilities.h"
#include "user/user_rate_functions.h"
#include "user/thermo.h"
#include "user/screen.h"
#include "user/nse_corr.h"

namespace po = boost::program_options;

#define S_USE_APPROXIMATE_WEAK_RATES "use approximate weak rates"

//##############################################################################
// get_input().
//##############################################################################

boost::tuple<Libnucnet *, std::string, std::string>
get_input( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  try
  {

    std::string s_nuc_xpath = "", s_reac_xpath = "", s_zone_xpath = "";

    std::string s_purpose = "\nPurpose: compute the entropy generation rate per nucleon for the input xml_file for the selected nuclei, reactions, and zones.";

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
       "zone_xpath",
       po::value<std::string>(),
       "XPath to select zones (default: all zones)"
      )
    ;

    po::variables_map vm;        
    po::store(po::parse_command_line( argc, argv, desc), vm );
    po::notify(vm);    

    if( argc < 2 || vm.count("help") == 1 )
    {
      std::cerr << "\nUsage: " << argv[0] << " xml_file [options]" << std::endl;
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

    if( vm.count("zone_xpath") == 1 )
    {
      s_zone_xpath = vm["zone_xpath"].as<std::string>();
    }

    p_nucnet =
      Libnucnet__new_from_xml( argv[1], "", "", s_zone_xpath.c_str() );

    if( !p_nucnet )
    {
      std::cerr << "Could not open file." << std::endl;
      exit( EXIT_FAILURE );
    }

    return
      boost::tuple<Libnucnet *, std::string, std::string>(
        p_nucnet, s_nuc_xpath, s_reac_xpath
      );

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

  boost::tuple<Libnucnet *, std::string, std::string> my_tuple;
  Libnucnet *p_my_nucnet;
  std::vector<nnt::Zone> zones;

  //============================================================================
  // Get input.
  //============================================================================

  my_tuple = get_input( argc, argv );

  p_my_nucnet = my_tuple.get<0>();

  //============================================================================
  // Get the zone list.
  //============================================================================
  
  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  //============================================================================
  // Update with approximate rates.
  //============================================================================

  if( (*(zone_list.begin())).hasProperty( S_USE_APPROXIMATE_WEAK_RATES ) )
  {
    user::aa522a25__update_net(
      Libnucnet__getNet( p_my_nucnet )
    );
  }

  //============================================================================
  // Register user rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Store zones.
  //============================================================================
  
  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
    zones.push_back( zone );
  }

  //============================================================================
  // Iterate the zones.
  //============================================================================
  
  for( size_t i = 0; i < zones.size(); i++ )
  {

    if( i > 0 )
    {
      Libnucnet__Zone__copy_net_views(
        zones[i].getNucnetZone(),
        zones[0].getNucnetZone()
      );
    }

    if( !zones[i].hasProperty( nnt::s_MU_NUE_KT ) )
    {
      zones[i].updateProperty(
        nnt::s_MU_NUE_KT,
        GSL_NEGINF
      );
    }

    user::update_rate_functions_data( zones[i] );

    if(
      zones[i].hasProperty( nnt::s_USE_SCREENING ) &&
      zones[i].getProperty<std::string>( nnt::s_USE_SCREENING ) == "yes"
    )
    {
      user::set_screening_function( zones[i] );
    }

    if(
      zones[i].hasProperty( nnt::s_USE_NSE_CORRECTION ) &&
      zones[i].getProperty<std::string>(
        nnt::s_USE_NSE_CORRECTION
      ) == "yes"
    )
    {
      user::set_nse_correction_function( zones[i] );
    }

    std::cout <<
      Libnucnet__Zone__getLabel( zones[i].getNucnetZone(), 1 ) << "  " <<
      Libnucnet__Zone__getLabel( zones[i].getNucnetZone(), 2 ) << "  " <<
      Libnucnet__Zone__getLabel( zones[i].getNucnetZone(), 3 ) << "  " <<
      user::compute_entropy_generation_rate(
        zones[i],
        zones[i].getNetView(
          (my_tuple.get<1>()).c_str(),
          (my_tuple.get<2>()).c_str()
        )
      ) << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

