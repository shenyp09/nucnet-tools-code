/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
//
//     This is free software; you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     Please see the src/README.txt file in this distribution for more
//     information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use nucnet routines to create
//       a new Libnucnet structure from an input xml file, 
//       supply user-defined screening and reverse ratio correction functions,
//       compute rates for a zone using those functions, toggle whether to
//       compute the reverse rate from detaile balance, print out the rates,
//       and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet.h>
#include "screen.h"
#include "coul_corr.h"
#include "user_rate_functions.h"

/*##############################################################################
// Prototypes.
//############################################################################*/

int
print_forward_and_reverse( Libnucnet__Reaction *, Libnucnet__Zone * );

int
print_forward_only( Libnucnet__Reaction *, Libnucnet__Zone * );

int
compare_reactions( const Libnucnet__Reaction *, const Libnucnet__Reaction * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;
  Libnucnet__NetView *p_view;

  /*============================================================================
  // User-supplied data.
  //==========================================================================*/

  double d_rhoe;
  struct user_coul_corr_data my_corr_data;
  struct user_screening_data my_screening_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc!= 8 ) {
      fprintf(
        stderr,
        "\nUsage: %s nucnet_file t9 rho label1 label2 label3 detailed\n\n",
        argv[0]
      );
      fprintf(
        stderr, "  nucnet_file = input nuclear network xml filename\n\n"
      );
      fprintf(
        stderr, "  t9 = input temperature (in 10^9 K)\n\n"
      );
      fprintf(
        stderr, "  rho = input mass density (in g/cc)\n\n"
      );
      fprintf(
        stderr, "  label1 = first zone label\n\n"
      );
      fprintf(
        stderr, "  label2 = second zone label\n\n"
      );
      fprintf(
        stderr, "  label3 = third zone label\n\n"
      );
      fprintf(
        stderr, "  detailed = set detailed balance (\"on\" or \"off\" )\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new_from_xml(
    argv[1], NULL, NULL, NULL
  );

  /*============================================================================
  // Register user-supplied rate functions.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Retrieve zone.
  //==========================================================================*/

  printf( "\n" );

  p_zone = Libnucnet__getZoneByLabels( 
    p_my_nucnet, argv[4], argv[5], argv[6] 
  );

  if( !p_zone ) {

    fprintf( stderr, "Zone not found!\n" );
    return EXIT_FAILURE;

  }

  /*============================================================================
  // Set screening data and function.
  //==========================================================================*/

  my_screening_data.dYe2 =
    Libnucnet__Zone__computeZMoment( p_zone, 2 );

  Libnucnet__Zone__setScreeningFunction(
    p_zone,
    (Libnucnet__Zone__screeningFunction) my_screening_function,
    &my_screening_data
  );

  /*============================================================================
  // Set correction function.
  //==========================================================================*/

  my_corr_data.dA1 = -0.9052;
  my_corr_data.dA2 =  0.6322;

  Libnucnet__Zone__setNseCorrectionFactorFunction(
    p_zone,
    (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
    &my_corr_data
  );

  /*============================================================================
  // Set user rate functions data.
  //==========================================================================*/

  d_rhoe =
    atof( argv[3] ) *
    Libnucnet__Zone__computeZMoment( p_zone, 1 );
  
  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    CF88_WEAK_FIT,
    &d_rhoe
  );

  /*============================================================================
  // Toggle reverse rate detailed balance.
  //==========================================================================*/

  Libnucnet__Zone__toggleReverseRateDetailedBalance( p_zone, argv[7] );

  /*============================================================================
  // Compute rates with screening and correction.
  //==========================================================================*/

  Libnucnet__Zone__computeRates( p_zone, atof( argv[2] ), atof( argv[3] ) );

  /*============================================================================
  // Print forward and reverse rates.
  //==========================================================================*/

  if( Libnucnet__Zone__isComputingReverseRatesFromDetailedBalance( p_zone ) )
  {
    printf( "\t\t\tReaction\t\t\t Forward     Reverse \n" );
    printf(
      "=======================================================\t ==========  ==========\n"
    );
  }
  else
  {
    printf( "\t\t\tReaction\t\t\t Forward\n" );
    printf(
      "=======================================================\t ==========\n"
    );
  }

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  p_view =
    Libnucnet__NetView__new(
      Libnucnet__Zone__getNet( p_zone ),
      NULL,
      NULL 
    );

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac(
      Libnucnet__NetView__getNet( p_view )
    ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  if( Libnucnet__Zone__isComputingReverseRatesFromDetailedBalance( p_zone ) )
  {
    Libnucnet__Reac__iterateReactions(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_view )
      ),
      (Libnucnet__Reaction__iterateFunction) print_forward_and_reverse,
      p_zone
    );
  }
  else
  {
    Libnucnet__Reac__iterateReactions(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_view )
      ),
      (Libnucnet__Reaction__iterateFunction) print_forward_only,
      p_zone
    );
  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Zone__clearScreeningFunction( p_zone );
  Libnucnet__Zone__clearNseCorrectionFactorFunction( p_zone );

  Libnucnet__NetView__free( p_view );

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// print_forward_and_reverse().
//############################################################################*/

int
print_forward_and_reverse(
  Libnucnet__Reaction *p_reaction, Libnucnet__Zone *p_zone
)
{

   double d_forward, d_reverse;

   Libnucnet__Zone__getRatesForReaction(
     p_zone,
     p_reaction,
     &d_forward,
     &d_reverse
   );

   printf( "%-55s%12.4e%12.4e\n",
     Libnucnet__Reaction__getString( p_reaction ),
     d_forward,
     d_reverse
   );

  return 1;

}

/*##############################################################################
// print_forward_only().
//############################################################################*/

int
print_forward_only(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Zone *p_zone
)
{

   double d_forward, d_reverse;

   Libnucnet__Zone__getRatesForReaction(
     p_zone,
     p_reaction,
     &d_forward,
     &d_reverse
   );

   printf( "%-55s%12.4e\n",
     Libnucnet__Reaction__getString( p_reaction ),
     d_forward
   );

  return 1;

}

/*##############################################################################
// compare_reactions()
//############################################################################*/

int
compare_reactions(
  const Libnucnet__Reaction *p_reaction1,
  const Libnucnet__Reaction *p_reaction2
)
{

  return
    strcmp(
      Libnucnet__Reaction__getString( p_reaction1 ),
      Libnucnet__Reaction__getString( p_reaction2 )
    );

}
