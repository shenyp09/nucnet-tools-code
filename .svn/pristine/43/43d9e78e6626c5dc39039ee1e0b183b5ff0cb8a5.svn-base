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
//! \file
//! \brief Example code to print weak forward and reverse rate for a
//!     reaction from input xml for given input.
////////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

#include <Libnucnet.h>
#include <Libstatmech.h>

#include "nnt/two_d_weak_rates.h"
#include "user/weak_utilities.h"

namespace po = boost::program_options;

//##############################################################################
// compute_rates().
//##############################################################################

std::pair<double, double>
compute_rates(
  Libnucnet__Net * p_net,
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  double d_rho,
  double d_ye,
  double d_mu_nue_kT
)
{
  
  double d_forward, d_reverse;

  double d_rhoe = d_ye * d_rho;

  Libstatmech__Fermion * p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  double d_muekT =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      d_t9 * GSL_CONST_NUM_GIGA,
      d_rhoe * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    );

  double d_eta_F =
    d_muekT
    +
    (
      Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( d_t9 )
    );

  if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      "two-d weak rates log10 ft"
    ) == 0
  )
  {
    d_forward =
      user::compute_weak_rate_from_log10_ft(
        p_reaction,
        p_net,
        Libstatmech__Fermion__getRestMass( p_electron ),
        d_t9,
        d_rhoe,
        d_eta_F,
        d_mu_nue_kT
      );
  }
  else if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      "two-d weak rates"
    ) == 0
  )
  {
    d_forward =
      user::compute_two_d_weak_rate(
        p_reaction,
        d_t9,
        &d_rhoe
      );
  }
  else
  {
    std::cerr << "Invalid weak reaction." << std::endl;
    exit( EXIT_FAILURE );
  }

  d_reverse =
    user::compute_reverse_weak_rate_for_reaction(
      Libnucnet__Net__getNuc( p_net ),
      p_reaction,
      d_forward,
      d_t9,
      d_rho,
      d_muekT,
      d_mu_nue_kT
    );

  Libstatmech__Fermion__free( p_electron );

  return std::make_pair( d_forward, d_reverse );

}

//##############################################################################
// at_option_parser().
//##############################################################################

std::pair<std::string, std::string>
at_option_parser(std::string const&s)
{
  if('@' == s[0])
    return std::make_pair(std::string("response-file"), s.substr(1));
  else
    return std::pair<std::string, std::string>();
}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet__Net * p_my_net;
  Libnucnet__Reaction * p_reaction;
  std::vector<double> t9_vector, rho_vector, ye_vector;
  double d_forward, d_reverse;
  double
    t9_min = 1., t9_max = 10., rho_min = 1.e7, rho_max = 1.e9,
    ye_min = 0.1, ye_max = 0.5, mu_nue_kT = GSL_NEGINF;
  size_t t9_points = 10, rho_points = 1, ye_points = 1;
  std::string s_t9_log = "linear", s_rho_log = "linear", s_ye_log = "linear";

  std::string s_purpose =
    "\nPurpose: compute two-d weak rates for reaction from reac_xml input file using input nuc_xml file.";

  try
  {

    po::options_description desc("\nOptions");
    desc.add_options()
      ( "help", "print out this help message and exit" )
      (
        "response-file",
        po::value<std::string>(),
        "can be specified with '@name', too"
      )
      (
       "t9_min",
       po::value<double>(),
       "minimum t9 (default: 1.)"
      )
      (
       "t9_max",
       po::value<double>(),
       "maximum t9 (default: 10.)"
      )
      (
       "t9_points",
       po::value<size_t>(),
       "number of t9 points (default: 10)"
      )
      (
       "t9_log",
       po::value<std::string>(),
       "linear or logarithmic spacing (default: linear)"
      )
      (
       "rho_min",
       po::value<double>(),
       "minimum rho (default: 1.e7)"
      )
      (
       "rho_max",
       po::value<double>(),
       "maximum rho (default: 1.e9)"
      )
      (
       "rho_points",
       po::value<size_t>(),
       "number of rho points (default: 1)"
      )
      (
       "rho_log",
       po::value<std::string>(),
       "linear or logarithmic spacing (default: linear)"
      )
      (
       "ye_min",
       po::value<double>(),
       "minimum electron-to-nucleon ratio (default: 0.1)"
      )
      (
       "ye_max",
       po::value<double>(),
       "maximum electron-to-nucleon ratio (default: 0.5)"
      )
      (
       "ye_points",
       po::value<size_t>(),
       "number of ye points (default: 1)"
      )
      (
       "ye_log",
       po::value<std::string>(),
       "linear or logarithmic spacing (default: linear)"
      )
      (
       "mu_nue_kT",
       po::value<double>(),
       "neutrino_e chemical potential / kT (default: -inf)"
      )
      ;

    po::variables_map vm;
    po::store( po::command_line_parser(argc, argv).options(desc)
              .extra_parser(at_option_parser).run(), vm );
    po::notify(vm);

    if( argc < 4 || vm.count("help") == 1 )
    {
      std::cerr <<
        "\nUsage: " <<
           argv[0] << " nuc_xml reac_xml reaction [options]" << std::endl;
      std::cerr << s_purpose << std::endl;
      std::cout << desc << "\n";
      exit( EXIT_FAILURE );
    }

    if( vm.count("response-file") )
    {
      std::ifstream ifs(vm["response-file"].as<std::string>().c_str());
      if( !ifs )
      {
        std::cout << "Could not open the response file\n";
        return EXIT_FAILURE;
      }
      std::stringstream ss;
      ss << ifs.rdbuf();
      boost::char_separator<char> sep(" \n\r");
      std::string sstr = ss.str();
      boost::tokenizer<boost::char_separator<char> > tok( sstr, sep );
      std::vector<std::string> args;
      std::copy( tok.begin(), tok.end(), back_inserter(args) );
      po::store( po::command_line_parser(args).options(desc).run(), vm );      
    }

    if( vm.count("t9_min") == 1 )
    {
      t9_min = vm["t9_min"].as<double>();
    }

    if( vm.count("t9_max") == 1 )
    {
      t9_max = vm["t9_max"].as<double>();
    }

    if( vm.count("rho_min") == 1 )
    {
      rho_min = vm["rho_min"].as<double>();
    }

    if( vm.count("rho_max") == 1 )
    {
      rho_max = vm["rho_max"].as<double>();
    }

    if( vm.count("ye_min") == 1 )
    {
      ye_min = vm["ye_min"].as<double>();
    }

    if( vm.count("ye_max") == 1 )
    {
      ye_max = vm["ye_max"].as<double>();
    }

    if( vm.count("t9_points") == 1 )
    {
      t9_points = vm["t9_points"].as<size_t>();
      if( t9_points == 0 )
      {
        std::cerr << "Number of t9 points must be > 0." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    if( vm.count("rho_points") == 1 )
    {
      rho_points = vm["rho_points"].as<size_t>();
      if( rho_points == 0 )
      {
        std::cerr << "Number of rho points must be > 0." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    if( vm.count("ye_points") == 1 )
    {
      ye_points = vm["ye_points"].as<size_t>();
      if( ye_points == 0 )
      {
        std::cerr << "Number of ye points must be > 0." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    if( vm.count("t9_log") == 1 )
    {
      s_t9_log = vm["t9_log"].as<std::string>();
      if( s_t9_log != "linear" && s_t9_log != "logarithmic" )
      {
        std::cerr <<
          "t9_log must be set as \"linear\" or \"logarithmic\"." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    if( vm.count("rho_log") == 1 )
    {
      s_rho_log = vm["rho_log"].as<std::string>();
      if( s_rho_log != "linear" && s_rho_log != "logarithmic" )
      {
        std::cerr <<
          "rho_log must be set as \"linear\" or \"logarithmic\"." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    if( vm.count("ye_log") == 1 )
    {
      s_ye_log = vm["ye_log"].as<std::string>();
      if( s_ye_log != "linear" && s_ye_log != "logarithmic" )
      {
        std::cerr <<
          "ye_log must be set as \"linear\" or \"logarithmic\"." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    if( vm.count("mu_nue_kT") == 1 )
    {
      mu_nue_kT = vm["mu_nue_kT"].as<double>();
    }

    t9_vector = nnt::get_vector( t9_min, t9_max, t9_points, s_t9_log );

    rho_vector = nnt::get_vector( rho_min, rho_max, rho_points, s_rho_log );

    ye_vector = nnt::get_vector( ye_min, ye_max, ye_points, s_ye_log );

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
  // Create network.
  //==========================================================================*/

  p_my_net = Libnucnet__Net__new();

  //============================================================================
  // Update nuclear and reaction data.
  //==========================================================================*/

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( p_my_net ), argv[1], NULL
  );

  Libnucnet__Reac__updateFromXml(
    Libnucnet__Net__getReac( p_my_net ), argv[2], NULL
  );

  //============================================================================
  // Get reaction.
  //==========================================================================*/

  p_reaction =
    Libnucnet__Reac__getReactionByString(
      Libnucnet__Net__getReac( p_my_net ),
      argv[3]
    );

  if( !p_reaction )
  {
    std::cerr << "Reaction " << argv[3] << " not found." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Loop over t9.
  //==========================================================================*/

  std::cout << std::endl << std::endl << "Reaction: " <<
    Libnucnet__Reaction__getString( p_reaction ) << std::endl << std::endl;

  boost::format fmt1(
    "    T9   \tRho(g/cc)\t    Ye   \t Forward\t Reverse\n"
  );
  std::cout << fmt1.str() << std::endl;

  boost::format
    fmt2( "==========\t==========\t==========\t==========\t==========\n" );
  std::cout << fmt2.str() << std::endl;

  for( size_t i_temp = 0; i_temp < t9_vector.size(); i_temp++ )
  {

    for( size_t i_rho = 0; i_rho < rho_vector.size(); i_rho++ )
    {

      for( size_t i_ye = 0; i_ye < ye_vector.size(); i_ye++ )
      {

        boost::tie( d_forward, d_reverse ) =
          compute_rates(
            p_my_net,
            p_reaction,
            t9_vector[i_temp],
            rho_vector[i_rho],
            ye_vector[i_ye],
            mu_nue_kT
          );

        fprintf(
          stdout,
          "%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
          t9_vector[i_temp],
          rho_vector[i_rho],
          ye_vector[i_ye],
          d_forward,
          d_reverse
        );

      }

    }
        
  }

  //============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Net__free( p_my_net );

  return EXIT_SUCCESS;

}

