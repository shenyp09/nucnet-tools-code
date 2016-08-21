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
//! \brief Example code to check whether a reaction changes sign of Q value
//!         with different nuclear mass data.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include <Libnucnet.h>
#include "nnt/auxiliary.h"

namespace po = boost::program_options;

//##############################################################################
// get_input().
//##############################################################################

std::pair<Libnucnet__Net *, Libnucnet__Net *>
get_input( int argc, char **argv )
{

  Libnucnet__Net * p_net_1, * p_net_2;

  try
  {

    std::string s_nuc_xpath = "", s_reac_xpath = "";

    std::string s_purpose = "\nPurpose: find reactions in reac_xml that have different signs for Q value for the nuclear data in nuc1_xml and nuc2_xml possibly selected by nuclide and reaction XPath.";

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

    if( argc < 2 || vm.count("help") == 1 )
    {
      std::cerr << "\nUsage: " << argv[0] << " nuc1_xml nuc2_xml reac_xml [options]" << std::endl;
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

    p_net_1 = Libnucnet__Net__new();

    p_net_2 = Libnucnet__Net__new();

    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( p_net_1 ),
      argv[1],
      s_nuc_xpath.c_str()
    );

    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( p_net_2 ),
      argv[2],
      s_nuc_xpath.c_str()
    );

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_net_1 ),
      argv[3],
      s_reac_xpath.c_str()
    );
  
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_net_2 ),
      argv[3],
      s_reac_xpath.c_str()
    );
  
    return std::make_pair( p_net_1, p_net_2 );

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

  Libnucnet__Net *p_net_1, *p_net_2;

  //============================================================================
  // Get input.
  //============================================================================

  boost::tie( p_net_1, p_net_2 ) = get_input( argc, argv );

  //============================================================================
  // Get reaction list.
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( p_net_1 ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( Libnucnet__Net__getReac( p_net_1 ) );

  //============================================================================
  // Print header.
  //============================================================================

  std::cout <<
    boost::format( "\n\t\t\tReaction\t\t\t\tFile 1  File 2\n" );
  std::cout <<
    boost::format(
      "=======================================================         ======  ======\n"
    );

  //============================================================================
  // Loop.
  //============================================================================

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    if(
      Libnucnet__Net__isValidReaction( p_net_1, reaction.getNucnetReaction() )
    )
    {
      double d_q1, d_q2;

      Libnucnet__Reaction * p_reaction =
	Libnucnet__Reac__getReactionByString(
	  Libnucnet__Net__getReac( p_net_2 ),
	  Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
	);

      d_q1 =
	Libnucnet__Net__computeReactionQValue(
	  p_net_1,
	  reaction.getNucnetReaction()
	);

      if( Libnucnet__Net__isValidReaction( p_net_2, p_reaction ) )
      {
        d_q2 =
          Libnucnet__Net__computeReactionQValue(
            p_net_2,
            p_reaction
        );

        if( d_q1 * d_q2 < 0 )
        {
          std::cout <<
            boost::format( "%-57s\t%6.3f\t%6.3f\n" ) %
            Libnucnet__Reaction__getString( reaction.getNucnetReaction() ) %
	    d_q1 %
	    d_q2;
        }
      }
      else
      {
        std::cout <<
          boost::format( "%-57s\t%6.3f\t------\n" ) %
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() ) %
	  d_q1;
      }
    }
     
  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Net__free( p_net_1 );
  Libnucnet__Net__free( p_net_2 );

  return EXIT_SUCCESS;

}
