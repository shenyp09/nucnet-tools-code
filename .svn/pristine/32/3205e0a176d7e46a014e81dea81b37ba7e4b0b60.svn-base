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
//! \brief Example code to convert text ffn-type rates file to xml file, storing
//!     both two-d weak rates and log10 ft values. 
////////////////////////////////////////////////////////////////////////////////
 
#include <boost/program_options.hpp>

#include <Libnucnet.h>
#include <Libstatmech.h>

#include "nnt/two_d_weak_rates.h"
#include "nnt/string_defs.h"
#include "nnt/param_defs.h"

#define MAX_SIZE 256
#define D_ALPHA    1.
#define D_MU_NUE_KT  GSL_NEGINF

namespace po = boost::program_options;

//##############################################################################
// Prototypes.
//##############################################################################

double
convert_to_log10_ft(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  Libstatmech__Fermion *,
  double,
  double,
  double,
  const char *
);

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet__Net *p_my_net;
  Libnucnet__Reaction *p_reaction_bp, *p_reaction_ec, *p_reaction_bm;
  Libnucnet__Reaction *p_reaction_pc;
  Libnucnet__Species *p_reactant, *p_target;
  Libstatmech__Fermion *p_electron;

  FILE *p_input_file;
  char s_source[MAX_SIZE], *s_end;
  int i_row, i_col, i_number_rows, i_number_columns;
  unsigned int i_z, i_a;
  char s_bp[MAX_SIZE], s_ec[MAX_SIZE], s_bm[MAX_SIZE], s_pc[MAX_SIZE];
  char s_t9[MAX_SIZE], s_log10_rhoe[MAX_SIZE];
  char s_row[MAX_SIZE], s_col[MAX_SIZE];
  bool b_convert = false;

  std::string s_purpose =
    "\nPurpose: convert text ffn-type rates text file reac_txt to xml file reac_xml.\n";

  try
  {

    po::options_description desc("\nOptions");
    desc.add_options()
      ( "help", "print out this help message and exit" )
      (
        "log10_ft",
        "convert rates to log10 ft where possible"
      )
      ;

    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );

    if( vm.count("help") || argc < 4 )
    {
      std::cout <<
        "\nUsage: " <<
           argv[0] << " nuc_xml reac_txt reac_xml [options]" << std::endl;
      std::cout << s_purpose << std::endl;
      std::cout << desc << "\n";
      return EXIT_FAILURE;
    }

    if( vm.count("log10_ft") )
    {
      b_convert = true;
    }

  }
  catch( std::exception& e )
  {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  catch(...)
  {
    std::cerr << "Exception of unknown type!\n";
  }

  //============================================================================
  // Create net.
  //============================================================================

  p_my_net = Libnucnet__Net__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( p_my_net ),
    argv[1],
    NULL
  );

  //============================================================================
  // Open input file.
  //============================================================================

  if( ( p_input_file = fopen( argv[2], "r" ) ) == NULL ) {
    fprintf( stderr, "Could not open file.\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read source and number of rows and columns.
  //============================================================================

  fgets( s_source, MAX_SIZE, p_input_file );

  s_end = strchr( s_source, '\n' );
  if( s_end ) {
    *s_end = '\0';
  }

  fscanf( p_input_file, "%d %d\n", &i_number_rows, &i_number_columns );

  //============================================================================
  // Get electron.
  //============================================================================

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  //============================================================================
  // Read in input.
  //============================================================================

  while( !feof( p_input_file ) ) {

  //----------------------------------------------------------------------------
  // Read in and get target.
  //----------------------------------------------------------------------------

    fscanf(
      p_input_file,
      "%u  %u\n",
      (unsigned int *) &i_z,
      (unsigned int *) &i_a
    );

    p_reactant =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z, i_a, NULL
      );

  //----------------------------------------------------------------------------
  // Create reactions and add reactants and products.
  //----------------------------------------------------------------------------

    //..........................................................................
    // beta plus.
    //..........................................................................

    p_reaction_bp = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_bp, s_source );

    Libnucnet__Reaction__addReactant(
      p_reaction_bp,
      Libnucnet__Species__getName( p_reactant )
    );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z - 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_bp,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_bp, "positron" );
    Libnucnet__Reaction__addProduct( p_reaction_bp, "neutrino_e" );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction_bp,
      nnt::s_TWO_D_WEAK_RATES
    );

    //..........................................................................
    // electron capture.
    //..........................................................................

    p_reaction_ec = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_ec, s_source );

    Libnucnet__Reaction__addReactant(
      p_reaction_ec,
      Libnucnet__Species__getName( p_reactant )
    );

    Libnucnet__Reaction__addReactant( p_reaction_ec, "electron" );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z - 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_ec,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_ec, "neutrino_e" );

    if( !b_convert )
    {
      Libnucnet__Reaction__setUserRateFunctionKey(
        p_reaction_ec,
        nnt::s_TWO_D_WEAK_RATES
      );
    }
    else
    {
      Libnucnet__Reaction__setUserRateFunctionKey(
        p_reaction_ec,
        nnt::s_TWO_D_WEAK_RATES_LOG10_FT
      );
    }

  //----------------------------------------------------------------------------
  // Reverse reactions. Decrease Z.
  //----------------------------------------------------------------------------

    p_reactant =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), --i_z, i_a, NULL
      );

    //..........................................................................
    // beta minus.
    //..........................................................................

    p_reaction_bm = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_bm, s_source );

    Libnucnet__Reaction__addReactant(
      p_reaction_bm,
      Libnucnet__Species__getName( p_reactant )
    );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z + 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_bm,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_bm, "electron" );
    Libnucnet__Reaction__addProduct( p_reaction_bm, "anti-neutrino_e" );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction_bm,
      nnt::s_TWO_D_WEAK_RATES
    );

    //..........................................................................
    // positron capture.
    //..........................................................................

    p_reaction_pc = Libnucnet__Reaction__new();

    Libnucnet__Reaction__updateSource( p_reaction_pc, s_source );

    Libnucnet__Reaction__addReactant(
      p_reaction_pc,
      Libnucnet__Species__getName( p_reactant )
    );

    Libnucnet__Reaction__addReactant( p_reaction_pc, "positron" );

    p_target =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_my_net ), i_z + 1, i_a, NULL
      );

    Libnucnet__Reaction__addProduct(
      p_reaction_pc,
      Libnucnet__Species__getName( p_target )
    );

    Libnucnet__Reaction__addProduct( p_reaction_pc, "anti-neutrino_e" );

    if( !b_convert )
    {
      Libnucnet__Reaction__setUserRateFunctionKey(
        p_reaction_pc,
        nnt::s_TWO_D_WEAK_RATES
      );
    }
    else
    {
      Libnucnet__Reaction__setUserRateFunctionKey(
        p_reaction_pc,
        nnt::s_TWO_D_WEAK_RATES_LOG10_FT
      );
    }


  //----------------------------------------------------------------------------
  // Read in data.
  //----------------------------------------------------------------------------

    for( i_col = 0; i_col < i_number_columns; i_col++ )
    {

      sprintf( s_col, "%d", i_col );

      for( i_row = 0; i_row < i_number_rows; i_row++ )
      {

        sprintf( s_row, "%d", i_row );

        fscanf(
          p_input_file,
          "%s %s %s %s %s %s\n",
          s_t9, s_log10_rhoe, s_bp, s_ec, s_bm, s_pc
        );  

    //..........................................................................
    // Update t9 tables for each reaction.
    //..........................................................................

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bp,
          nnt::s_T9,
          s_row,
          NULL,
          s_t9
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_T9,
          s_row,
          NULL,
          s_t9
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bm,
          nnt::s_T9,
          s_row,
          NULL,
          s_t9
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_T9,
          s_row,
          NULL,
          s_t9
        ); 

    //..........................................................................
    // Update rhoe tables for each reaction.
    //..........................................................................

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bp,
          nnt::s_LOG10_RHOE,
          s_col,
          NULL,
          s_log10_rhoe
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_LOG10_RHOE,
          s_col,
          NULL,
          s_log10_rhoe
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bm,
          nnt::s_LOG10_RHOE,
          s_col,
          NULL,
          s_log10_rhoe
        ); 

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_LOG10_RHOE,
          s_col,
          NULL,
          s_log10_rhoe
        ); 

    //..........................................................................
    // Update rate tables for each reaction.
    //..........................................................................

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bp,
          nnt::s_LOG10_RATE,
          s_row,
          s_col,
          s_bp
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_ec,
          nnt::s_LOG10_RATE,
          s_row,
          s_col,
          s_ec
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_bm,
          nnt::s_LOG10_RATE,
          s_row,
          s_col,
          s_bm
        );

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          p_reaction_pc,
          nnt::s_LOG10_RATE,
          s_row,
          s_col,
          s_pc
        );

    //..........................................................................
    // Add log10 ft to reactions.
    //..........................................................................

        if( b_convert )
        {

	  Libnucnet__Reaction__updateUserRateFunctionProperty(
	    p_reaction_ec,
	    nnt::s_LOG10_FT,
	    s_row,
	    s_col,
	    boost::lexical_cast<std::string>(
	      convert_to_log10_ft(
		p_reaction_ec,
		p_my_net,
		p_electron,
		atof( s_t9 ),
		pow( 10., atof( s_log10_rhoe ) ), 
		D_MU_NUE_KT,
		s_ec
	      )
	    ).c_str()
	  );

	  Libnucnet__Reaction__updateUserRateFunctionProperty(
	    p_reaction_pc,
	    nnt::s_LOG10_FT,
	    s_row,
	    s_col,
	    boost::lexical_cast<std::string>(
	      convert_to_log10_ft(
		p_reaction_pc,
		p_my_net,
		p_electron,
		atof( s_t9 ),
		pow( 10., atof( s_log10_rhoe ) ), 
		D_MU_NUE_KT,
		s_pc
	      )
	    ).c_str()
	  );

        }

      }

    }

  //----------------------------------------------------------------------------
  // Add reactions.
  //----------------------------------------------------------------------------

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_bp
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_ec
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_bm
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_my_net ),
      p_reaction_pc
    );

  }     
  
  //============================================================================
  // Close input file.
  //============================================================================

  fclose( p_input_file );

  //============================================================================
  // Write to XML.
  //============================================================================

  Libnucnet__Reac__writeToXmlFile(
    Libnucnet__Net__getReac( p_my_net ),
    argv[3]
  );

  //============================================================================
  // Clean up.
  //============================================================================

  Libstatmech__Fermion__free( p_electron );
  Libnucnet__Net__free( p_my_net );

  //============================================================================
  // Done!
  //============================================================================

  return EXIT_SUCCESS;

}

//##############################################################################
// convert_to_log10_ft().
//##############################################################################

double
convert_to_log10_ft(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Net *p_net,
  Libstatmech__Fermion *p_electron,
  double d_t9,
  double d_rhoe,
  double d_mu_nue_kT,
  const char *s_value
)
{

  double d_I, d_eta_F;

  d_eta_F = 
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      d_t9 * GSL_CONST_NUM_GIGA,
      d_rhoe * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    ) +
    Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( d_t9 );

  d_I = 
    nnt::ffnIV__compute_Ie(
      p_reaction,
      p_net, 
      Libstatmech__Fermion__getRestMass( p_electron ),
      d_t9,
      d_eta_F,
      d_mu_nue_kT
    );

  if( d_I > 0. )
    return
      log10(
        M_LN2
      ) +
      log10( d_I ) -
      atof( s_value ) +
      log10( D_ALPHA );
  else
    return 0;

}
