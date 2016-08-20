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
//       Example to demonstrate how to use libnucnet routines to
//       validate an input Libnucnet__Reac xml file against Webnucleo's
//       schema over the web.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>

int main( int argc, char * argv[] ) {

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc != 2 ) {
      fprintf(
        stderr, "\nUsage: %s reac_filename\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_filename = input reactions xml filename\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Validate reaction file
  //==========================================================================*/

  if( Libnucnet__Reac__is_valid_input_xml( argv[1] ) )
    fprintf( stdout, "Valid input Libnucnet__Reac xml file!\n" );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

