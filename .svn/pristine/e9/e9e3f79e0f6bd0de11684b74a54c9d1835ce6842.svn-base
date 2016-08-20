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
//!
//! \file
//! \brief Example code to compute the energy generation rate 
//!    in zones in a output network xml file.
//!
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include <numeric>

#include <boost/program_options.hpp>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

#include "user/user_rate_functions.h"
#include "user/flow_utilities.h"
#include "user/aa522a25.h"
#include "user/thermo.h"

#define S_USE_APPROXIMATE_WEAK_RATES "use approximate weak rates"

namespace po = boost::program_options;

//##############################################################################
// Types.
//##############################################################################

typedef std::vector<std::pair<std::string, double> > result_t;

//##############################################################################
// my_sorter()>
//##############################################################################

bool
my_sorter(
  std::pair<std::string, double> p1,
  std::pair<std::string, double> p2
)
{
  return fabs( p1.second ) > fabs( p2.second );
}

//##############################################################################
// get_input().
//##############################################################################

boost::tuple<Libnucnet *, std::string, std::string, double>
get_input( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  try
  {

    std::string s_nuc_xpath = "", s_reac_xpath = "", s_zone_xpath = "";
    double d_cutoff = 1.e-10;

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
      (
       "cutoff",
       po::value<double>(),
       "cutoff for reaction printout (default: 1.e-10)"
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

    if( vm.count("cutoff") == 1 )
    {
      d_cutoff = vm["cutoff"].as<double>();
    }

    p_nucnet =
      Libnucnet__new_from_xml( argv[1], "", "", s_zone_xpath.c_str() );

    if( !p_nucnet )
    {
      std::cerr << "Could not open file." << std::endl;
      exit( EXIT_FAILURE );
    }

    return
      boost::tuple<Libnucnet *, std::string, std::string, double>(
        p_nucnet, s_nuc_xpath, s_reac_xpath, d_cutoff
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

  std::vector<nnt::Zone> zones;
  Libnucnet__NetView * p_view;
  std::string s_nuc_xpath, s_reac_xpath;
  boost::tuple<Libnucnet *, std::string, std::string, double> my_tuple;

  //============================================================================
  // Check input.
  //============================================================================

  my_tuple = get_input( argc, argv );

  //============================================================================
  // Get a view.
  //============================================================================
  
  p_view =
    Libnucnet__NetView__new(
      Libnucnet__getNet( my_tuple.get<0>() ),
      my_tuple.get<1>().c_str(),
      my_tuple.get<2>().c_str()
    );

  //============================================================================
  // Get the zones.
  //============================================================================
  
  Libnucnet__setZoneCompareFunction(
    my_tuple.get<0>(),
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( my_tuple.get<0>() );

  //============================================================================
  // Update with approximate rates.
  //============================================================================

  if( (*(zone_list.begin())).hasProperty( S_USE_APPROXIMATE_WEAK_RATES ))
  {
    user::aa522a25__update_net(
      Libnucnet__getNet( my_tuple.get<0>() )
    );
  }

  //============================================================================
  // Register user rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( my_tuple.get<0>() ) )
  );

  //============================================================================
  // Iterate the zones.
  //============================================================================

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    zones.push_back( zone );

  }

  std::vector<result_t> results( zones.size() );

  for( size_t i = 0; i < zones.size(); i++ )
  {

    results[i] =
      user::compute_zone_reactions_entropy_generation_rate_per_nucleon(
        zones[i],
        p_view
    );

    std::sort( results[i].begin(), results[i].end(), my_sorter );

  }

  //============================================================================
  // Output.
  //============================================================================

  for( size_t i = 0; i < zones.size(); i++ )
  {

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

    std::vector<double> v_sdot;
    for( size_t j = 0; j < results[i].size(); j++ )
    {
      v_sdot.push_back( (results[i])[j].second );
    }

    std::cout <<
      boost::format( "time(s) = %g  t9 = %g rho(g/cc) = %g sdot = %g\n\n" ) %
        zones[i].getProperty<double>( nnt::s_TIME ) %
        zones[i].getProperty<double>( nnt::s_T9 ) %
        zones[i].getProperty<double>( nnt::s_RHO ) %
        std::accumulate( v_sdot.begin(), v_sdot.end(), 0. );

    for( size_t j = 0; j < results[i].size(); j++ )
    {
      if( fabs( (results[i])[j].second ) >= my_tuple.get<3>() )
      {
        std::cout << 
          boost::format( "%50s  %e\n" ) %
            (results[i])[j].first %
            (results[i])[j].second;
      }
    }

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__NetView__free( p_view );
  Libnucnet__free( my_tuple.get<0>() );

  return EXIT_SUCCESS;

}
