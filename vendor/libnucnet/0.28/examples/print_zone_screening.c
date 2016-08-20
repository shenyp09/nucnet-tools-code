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
//       compute rates for a zone using those functions, print out
//       the screening and correction factors, and clear the
//       structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet.h>
#include "screen.h"
#include "coul_corr.h"
#include "user_rate_functions.h"

#define  S_T9      "t9"
#define  S_RHO     "rho"
#define  S_YE      "Ye"

/*##############################################################################
// Prototypes.
//############################################################################*/

int
print_reaction( Libnucnet__Reaction *, Libnucnet__Zone * );

int
compare_reactions( const Libnucnet__Reaction *, const Libnucnet__Reaction * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;
  Libnucnet__NetView *p_view;
  char s_ye[32];
  double d_ye;

  /*============================================================================
  // User-supplied data.
  //==========================================================================*/

  struct user_coul_corr_data my_corr_data;
  struct user_screening_data my_screening_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc!= 7 ) {
      fprintf(
        stderr,
        "\nUsage: %s nucnet_file t9 rho label1 label2 label3\n\n",
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
  // Set data.
  //==========================================================================*/

  Libnucnet__Zone__updateProperty( p_zone, S_T9, NULL, NULL, argv[2] );

  Libnucnet__Zone__updateProperty( p_zone, S_RHO, NULL, NULL, argv[3] );

  d_ye = Libnucnet__Zone__computeZMoment( p_zone, 1 );

  sprintf( s_ye, "%f", d_ye );

  Libnucnet__Zone__updateProperty( p_zone, S_YE, NULL, NULL, s_ye );

  /*============================================================================
  // Print screening and correction factors.
  //==========================================================================*/

  printf( "\t\t\tReaction\t\t\t For. Screen  Rev. Screen Correction \n" );
  printf(
    "=======================================================\t ===========  =========== ==========\n"
  );

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

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac(
      Libnucnet__NetView__getNet( p_view )
    ),
    (Libnucnet__Reaction__iterateFunction) print_reaction,
    p_zone
  );

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
// print_reaction().
//############################################################################*/

int
print_reaction( Libnucnet__Reaction *p_reaction, Libnucnet__Zone *p_zone )
{

   double d_screen_f, d_screen_r, d_nse_corr;

   double d_t9, d_rho, d_Ye;

   d_t9 =
     atof(
       Libnucnet__Zone__getProperty( p_zone, S_T9, NULL, NULL )
     );

   d_rho =
     atof(
       Libnucnet__Zone__getProperty( p_zone, S_RHO, NULL, NULL )
     );

   d_Ye =
     atof(
       Libnucnet__Zone__getProperty( p_zone, S_YE, NULL, NULL )
     );

   my_reaction_screening_function(
     Libnucnet__Zone__getNet( p_zone ),
     p_reaction,
     d_t9,
     d_rho,
     d_Ye,
     Libnucnet__Zone__getScreeningData( p_zone ),
     &d_screen_f,
     &d_screen_r
   );
 
   d_nse_corr =
     Libnucnet__Net__computeReverseRatioCorrectionFactorForReaction(
       Libnucnet__Zone__getNet( p_zone ),
       p_reaction,
       d_t9,
       d_rho,
       d_Ye,
       Libnucnet__Zone__getNseCorrectionFactorFunction( p_zone ),
       Libnucnet__Zone__getNseCorrectionFactorData( p_zone )
     );

   printf( "%-55s%13.4e%13.4e%11.4e\n",
     Libnucnet__Reaction__getString( p_reaction ),
     d_screen_f,
     d_screen_r,
     d_nse_corr
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
