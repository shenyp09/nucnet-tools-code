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
//! \brief Example code to compare various equilibria to a network calculation.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#ifndef NO_OPENMP
#include <omp.h>
#endif

#include <Libnucnet.h>
#include <Libnuceq.h>

#include <boost/program_options.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/math.h"
#include "nnt/string_defs.h"
#include "user/aa522a25.h"
#include "user/nse_corr.h"

#include "user/weak_utilities.h"
#include "user/network_utilities.h"

#define S_WSE     "wse"
#define S_XPATH   "[z > 2]"

namespace po = boost::program_options;

//##############################################################################
// get_input().
//##############################################################################

boost::tuple<Libnucnet *, std::vector<std::string>, std::string>
get_input( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  try
  {

    std::string
      s_zone_xpath = "", s_cluster_xpath = S_XPATH, s_compute_wse = "";
    std::vector<std::string> s_clusters;

    std::string s_purpose = "\nPurpose: compare the network abundances in a network xml_file for the zones and equilibria.";

    std::string s_cluster_desc =
      "XPath to select cluster(s) for QSE calculation (default: \"" + std::string( S_XPATH ) + "\")";

    po::options_description desc("\nAllowed options");
    desc.add_options()
      ( "help", "print out this help message and exit" )
      (
       "zone_xpath",
       po::value<std::string>(),
       "XPath to select zones (default: all zones)"
      )
      (
       "cluster",
       po::value< std::vector<std::string> >()->composing(),
       s_cluster_desc.c_str()
      )
      (
       "wse",
       "also compute weak statistical equilibrium (default: do not compute)"
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

    if( vm.count("zone_xpath") == 1 )
    {
      s_zone_xpath = vm["zone_xpath"].as<std::string>();
    }

    if( vm.count("cluster") )
    {
      s_clusters = vm["cluster"].as<std::vector< std::string> >();
    }
    else
    {
      s_clusters.push_back( S_XPATH );
    }

    if( vm.count("wse") == 1 )
    {
      s_compute_wse = S_WSE;
    }

    p_nucnet =
      Libnucnet__new_from_xml( argv[1], "", "", s_zone_xpath.c_str() );

    if( !p_nucnet )
    {
      std::cerr << "Could not open file." << std::endl;
      exit( EXIT_FAILURE );
    }

    return
      boost::tuple<Libnucnet *, std::vector<std::string>, std::string>(
        p_nucnet, s_clusters, s_compute_wse
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
// coul_wse_function().
//##############################################################################

double
coul_wse_function(
  Libnuceq * p_equil,
  user::nse_corr_data_t * p_data
)
{

  double d_result =
    user::nse_correction(
      Libnucnet__Nuc__getSpeciesByName(
        Libnuceq__getNuc( p_equil ),
        "h1"
      ),
      Libnuceq__getT9( p_equil ),
      Libnuceq__getRho( p_equil ),
      Libnuceq__computeZMoment( p_equil, 1 ),
      p_data
    );

  return d_result;

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  std::vector<nnt::Zone> zones;
  std::vector<std::string> s_clusters;
  std::string s_compute_wse;
  std::vector<double> yh, yh_nse, yh_wse, ye_wse;

  //============================================================================
  // Get input.
  //============================================================================

  boost::tie( p_my_nucnet, s_clusters, s_compute_wse ) =
    get_input( argc, argv );

  //============================================================================
  // Check zones.
  //============================================================================

  if( Libnucnet__getNumberOfZones( p_my_nucnet ) == 0 )
  {
    std::cerr << "No zones." << std::endl;
    return EXIT_FAILURE;
  }
 
  //============================================================================
  // Output largest Z.
  //============================================================================

  fprintf(
    stdout,
    "Largest Z = %d\n",
    Libnucnet__Nuc__getLargestNucleonNumber(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      "z"
    )
  );

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

    zones.push_back( zone );
    yh_nse.push_back( 0. );
    yh_wse.push_back( 0. );
    ye_wse.push_back( 0. );
    yh.push_back( 0. );

    if( argc == 4 )
    {
      user::register_rate_functions(
        Libnucnet__Net__getReac(
          Libnucnet__Zone__getNet(
            zone.getNucnetZone()
          )
        )
      );

      user::update_rate_functions_data( zone );
    } 

  }

  zones[0].getNetView( S_XPATH, "" );

  for( size_t i = 1; i < zones.size(); i++ )
  {
    Libnucnet__Zone__copy_net_views(
      zones[i].getNucnetZone(),
      zones[0].getNucnetZone()
    );
  }

  user::set_weak_views_in_zones( p_my_nucnet );

  std::vector<gsl_vector *> z_net( zones.size() );
  std::vector<gsl_vector *> z_nse( zones.size() );
  std::vector<gsl_vector *> z_qse( zones.size() );
  std::vector<gsl_vector *> z_wse( zones.size() );

  #pragma omp parallel for schedule( dynamic, 1 )
  for( size_t i = 0; i < zones.size(); i++ )
  {
 
    boost::any my_nse_corr_data;

    gsl_vector * p_net_abundances =
      Libnucnet__Zone__getAbundances( zones[i].getNucnetZone() );

    z_net[i] =
      Libnucnet__Zone__getSummedAbundances( zones[i].getNucnetZone(), "z" );

    yh[i] =
      user::compute_cluster_abundance_moment(
        zones[i],
        S_XPATH,
        "a",
        0
      );

    //==========================================================================
    // Compute NSE.
    //==========================================================================

    Libnuceq * p_equil =
      Libnuceq__new(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
      );

    //==========================================================================
    // Add Coulomb correction if present.
    //==========================================================================
    
    if(
      zones[i].hasProperty( nnt::s_USE_NSE_CORRECTION )
      &&
      zones[i].getProperty<std::string>( nnt::s_USE_NSE_CORRECTION) == "yes"
    )
    {

      user::set_nse_correction_function( zones[i] );

      my_nse_corr_data =
        boost::any_cast<boost::function<boost::any()> >(
          zones[i].getFunction( nnt::s_NSE_CORRECTION_FACTOR_DATA_FUNCTION )
        )();

      Libnuceq__setNseCorrectionFactorFunction(
        p_equil,
        (Libnucnet__Species__nseCorrectionFactorFunction)
          user::nse_correction,
        &my_nse_corr_data
      );

    }

    //==========================================================================
    // Compute NSE.
    //==========================================================================

    Libnuceq__setYe(
      p_equil,
      Libnucnet__Zone__computeZMoment( zones[i].getNucnetZone(), 1 )
    );

    Libnuceq__computeEquilibrium(
      p_equil, 
      zones[i].getProperty<double>( nnt::s_T9 ),
      zones[i].getProperty<double>( nnt::s_RHO )
    );

    gsl_vector * p_abundances = Libnuceq__getAbundances( p_equil );

    Libnucnet__Zone__updateAbundances(
      zones[i].getNucnetZone(),
      p_abundances
    );

    z_nse[i] =
      Libnucnet__Zone__getSummedAbundances( zones[i].getNucnetZone(), "z" );

    yh_nse[i] =
      user::compute_cluster_abundance_moment(
        zones[i],
        S_XPATH,
        "a",
        0
      );

    gsl_vector_free( p_abundances );

    Libnucnet__Zone__updateAbundances(
      zones[i].getNucnetZone(),
      p_net_abundances
    );

    //==========================================================================
    // Compute QSE.
    //==========================================================================

    for( size_t j = 0; j < s_clusters.size(); j++ )
    {

      Libnuceq__Cluster__updateConstraint(
        Libnuceq__newCluster( p_equil, s_clusters[j].c_str() ),
        user::compute_cluster_abundance_moment(
          zones[i],
          s_clusters[j].c_str(),
          "a",
          0
        )
      );

    }

    Libnuceq__computeEquilibrium(
      p_equil, 
      zones[i].getProperty<double>( nnt::s_T9 ),
      zones[i].getProperty<double>( nnt::s_RHO )
    );

    p_abundances = Libnuceq__getAbundances( p_equil );

    Libnucnet__Zone__updateAbundances(
      zones[i].getNucnetZone(),
      p_abundances
    );

    z_qse[i] = 
      Libnucnet__Zone__getSummedAbundances(
        zones[i].getNucnetZone(),
        "z"
      );

    gsl_vector_free( p_abundances );

    for( size_t j = 0; j < s_clusters.size(); j++ )
    {
      Libnuceq__removeCluster(
        p_equil,
        Libnuceq__getCluster( p_equil, s_clusters[j].c_str() )
      );
    }

    Libnucnet__Zone__updateAbundances(
      zones[i].getNucnetZone(),
      p_net_abundances
    );

    //==========================================================================
    // Compute WSE.  Do it last for a zone because it modifies Ye.
    //==========================================================================

    if( s_compute_wse == S_WSE )
    {

      if(
	gsl_isinf(
          user::compute_thermo_quantity(
            zones[i],
            nnt::s_CHEMICAL_POTENTIAL_KT,
            nnt::s_NEUTRINO_E
          )
	) != -1
      )
      {

	Libnuceq__clearYe( p_equil );

	Libnuceq__computeEquilibrium(
	  p_equil, 
	  zones[i].getProperty<double>( nnt::s_T9 ),
	  zones[i].getProperty<double>( nnt::s_RHO )
	);

	p_abundances = Libnuceq__getAbundances( p_equil );

	ye_wse[i] = Libnuceq__computeZMoment( p_equil, 1 );

	Libnucnet__Zone__updateAbundances(
          zones[i].getNucnetZone(), p_abundances
        );

	yh_wse[i] = 
	  user::compute_cluster_abundance_moment(
	    zones[i],
	    S_XPATH,
	    "a",
	    0
	  );

	z_wse[i] =
	  Libnucnet__Zone__getSummedAbundances( zones[i].getNucnetZone(), "z" );

	gsl_vector_free( p_abundances );

      }
      else
      {

	ye_wse[i] =
	  nnt::compute_1d_root(
            boost::bind( user::yedot_root, _1, boost::ref( zones[i] ) ),
            0.5,
            1.1
          );

	yh_wse[i] =
	  user::compute_cluster_abundance_moment(
	    zones[i],
	    S_XPATH,
	    "a",
	    0
	  );

	z_wse[i] =
	  Libnucnet__Zone__getSummedAbundances( zones[i].getNucnetZone(), "z" );

      }

    }

    Libnucnet__Zone__updateAbundances(
      zones[i].getNucnetZone(),
      p_net_abundances
    );

    gsl_vector_free( p_net_abundances );
    Libnucnet__Zone__clearNetViews( zones[i].getNucnetZone() );
    Libnucnet__Zone__clearRates( zones[i].getNucnetZone() );

    Libnuceq__free( p_equil );

  }


  for( size_t i = 0; i < zones.size(); i++ )
  {

    //==========================================================================
    // Print out.
    //==========================================================================

    if( zones[i].hasProperty( nnt::s_TIME ) )
      std::cout << "time = " << zones[i].getProperty<double>( nnt::s_TIME ) <<
      std::endl;

    std::cout << "t9 = " << zones[i].getProperty<double>( nnt::s_T9 ) <<
      std::endl;

    std::cout << "rho = " << zones[i].getProperty<double>( nnt::s_RHO ) <<
      std::endl;

    std::cout << "Ye = " << 
      Libnucnet__Zone__computeZMoment( zones[i].getNucnetZone(), 1 ) <<
      std::endl;

    if( s_compute_wse == S_WSE )
      std::cout << "Ye_wse = " << ye_wse[i] << std::endl;

    std::cout << "Yh = " << yh[i] << std::endl;

    std::cout << "Yh_nse = " << yh_nse[i] << std::endl;

    if( s_compute_wse == S_WSE )
      std::cout << "Yh_wse = " << yh_wse[i] << std::endl;

    for(
      size_t j = 0;
      j <= Libnucnet__Nuc__getLargestNucleonNumber(
	     Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
	     "z"
	   );
      j++
    )
    {

      if( s_compute_wse == S_WSE )
      {
	fprintf(
	  stdout,
	  "%lu %e %e %e %e\n",
	  (unsigned long) j,
	  gsl_vector_get( z_net[i], j ),
	  gsl_vector_get( z_nse[i], j ),
	  gsl_vector_get( z_qse[i], j ),
	  gsl_vector_get( z_wse[i], j )
	);
      }
      else
      {
	fprintf(
	  stdout,
	  "%lu %e %e %e\n",
	  (unsigned long) j,
	  gsl_vector_get( z_net[i], j ),
	  gsl_vector_get( z_nse[i], j ),
	  gsl_vector_get( z_qse[i], j )
	);
      }

    }

    gsl_vector_free( z_net[i] );
    gsl_vector_free( z_nse[i] );
    gsl_vector_free( z_qse[i] );
    gsl_vector_free( z_wse[i] );

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
