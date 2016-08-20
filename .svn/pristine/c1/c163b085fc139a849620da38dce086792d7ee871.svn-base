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
//! \brief Example code to compute the screened and coulomb-corrected rates
//!    in zones in a output network xml file.
//!
////////////////////////////////////////////////////////////////////////////////

#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

#include <Libnucnet.h>

#include "nnt/wrappers.hpp"
#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

#include "user/aa522a25.h"
#include "user/user_rate_functions.h"
#include "user/screen.h"
#include "user/nse_corr.h"

namespace po = boost::program_options;

//##############################################################################
// get_input().
//##############################################################################

std::pair<Libnucnet *, double>
get_input( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  try
  {

    std::string s_nuc_xpath = "", s_reac_xpath = "", s_zone_xpath = "";
    double d_cutoff = 1.e-25;

    std::string s_purpose = "\nPurpose: compute the screened and unscreened reaction rates for zones in xml_file for the selected nuclei, reactions, and zones.";

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
       "minimum rate value to print out (default: 1.e-25)"
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
      if( d_cutoff <= 0 )
      {
        std::cerr << "Rate cutoff must be > 0." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    p_nucnet =
      Libnucnet__new_from_xml(
        argv[1],
        s_nuc_xpath.c_str(),
        s_reac_xpath.c_str(),
        s_zone_xpath.c_str()
      );

    if( !p_nucnet )
    {
      std::cerr << "Could not open file." << std::endl;
      exit( EXIT_FAILURE );
    }

    return std::make_pair( p_nucnet, d_cutoff );

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

  std::pair<Libnucnet *, double> my_pair;
  std::map<std::string, std::pair<double,double> > my_map;
  double d_forward, d_reverse;

  //============================================================================
  // Get input.
  //============================================================================

  my_pair = get_input( argc, argv );

  //============================================================================
  // Get the zone list.
  //============================================================================
  
  Libnucnet__setZoneCompareFunction(
    my_pair.first,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( my_pair.first );

  //============================================================================
  // Update with approximate rates.
  //============================================================================

  if(
    (*(zone_list.begin())).hasProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) &&
    (*(zone_list.begin())).getProperty<std::string>(
      nnt::s_USE_APPROXIMATE_WEAK_RATES
    ) == "yes"
  )
  {
    user::aa522a25__update_net(
      Libnucnet__getNet( my_pair.first )
    );
  }

  //============================================================================
  // Register user rate functions.
  //============================================================================

  user::register_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( my_pair.first ) )
  );

  //============================================================================
  // Order the reactions and get the reaction list.
  //============================================================================
  
  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( Libnucnet__getNet( my_pair.first ) ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__getNet( my_pair.first ) )
    );

  //============================================================================
  // Get printout formats.
  //============================================================================

  boost::format fmt1(
    "time(s) = %9.4e  t9 = %6.4f  rho(g/cc) = %9.4e  Ye = %6.4f"
  );

  boost::format fmt2(
    "t9 = %6.4f  rho(g/cc) = %9.4e  Ye = %6.4f"
  );

  boost::format fmt3(
    "\n\t\t\t        \t\t\t      Screened Rates                Unscreened Rates\n"
  );

  boost::format fmt4(
    "\n\t\t\tReaction\t\t\t  Forward Rate   Reverse Rate  Forward Rate   Reverse Rate\n"
  );

  boost::format fmt5(
   "=======================================================   ============   ============  ============   ============\n"
  );

  boost::format fmt6( "%-55s%14.5e%15.5e%14.5e%15.5e" );

  //============================================================================
  // Iterate the zones.
  //============================================================================
  
  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    boost::any screening_data, nse_corr_data;

    //==========================================================================
    // Print conditions.
    //==========================================================================

    if( zone.hasProperty( nnt::s_TIME ) )
    {
      fmt1 %
        zone.getProperty<double>( nnt::s_TIME ) %
        zone.getProperty<double>( nnt::s_T9 ) %
        zone.getProperty<double>( nnt::s_RHO ) %
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );
      std::cout << fmt1.str() << std::endl;
    }
    else
    {
      fmt2 %
        zone.getProperty<double>( nnt::s_TIME ) %
        zone.getProperty<double>( nnt::s_T9 ) %
        zone.getProperty<double>( nnt::s_RHO );
      std::cout << fmt2.str() << std::endl;
    }

    std::cout << fmt3.str() << std::endl;

    std::cout << fmt4.str() << std::endl;

    std::cout << fmt5.str() << std::endl;

    if( !zone.hasProperty( nnt::s_MU_NUE_KT ) )
    {
      zone.updateProperty(
        nnt::s_MU_NUE_KT,
        GSL_NEGINF
      );
    }

    user::update_rate_functions_data( zone );

    // Compute the unscreened rates.

    Libnucnet__Zone__computeRates(
      zone.getNucnetZone(),
      zone.getProperty<double>( nnt::s_T9 ),
      zone.getProperty<double>( nnt::s_RHO )
    );

    BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
    {

      Libnucnet__Zone__getRatesForReaction(
        zone.getNucnetZone(),
        reaction.getNucnetReaction(),
        &d_forward,
        &d_reverse
      );    

      my_map[Libnucnet__Reaction__getString( reaction.getNucnetReaction() )] =
        std::make_pair( d_forward, d_reverse );

    }

    // Compute the screened rates.
    
    user::set_screening_function( zone );

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

    user::set_nse_correction_function( zone );

    nse_corr_data =
      boost::any_cast<boost::function<boost::any()> >(
        zone.getFunction( nnt::s_NSE_CORRECTION_FACTOR_DATA_FUNCTION )
      )();

    Libnucnet__Zone__setNseCorrectionFactorFunction(
      zone.getNucnetZone(),
      (Libnucnet__Species__nseCorrectionFactorFunction)
        Libnucnet__Zone__getNseCorrectionFactorFunction( zone.getNucnetZone() ),
      &nse_corr_data
    );

    Libnucnet__Zone__computeRates(
      zone.getNucnetZone(),
      zone.getProperty<double>( nnt::s_T9 ),
      zone.getProperty<double>( nnt::s_RHO )
    );

    BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
    {

      Libnucnet__Zone__getRatesForReaction(
        zone.getNucnetZone(),
        reaction.getNucnetReaction(),
        &d_forward,
        &d_reverse
      );    
        
      if( d_forward > my_pair.second || d_reverse > my_pair.second )
      {

        
        std::string s_reac =
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() );

        fmt6 %
          s_reac %
          d_forward %
          d_reverse %
          my_map[s_reac].first %
          my_map[s_reac].second;

        std::cout << fmt6.str() << std::endl;

      }

    }

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( my_pair.first );

  return EXIT_SUCCESS;

}
