/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by George C. Jordan, IV and
//     Bradley S. Meyer.
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
//       Example of a user-supplied Coulomb correction factor function
//       header file.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet.h>

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// User-supplied data structure.
//############################################################################*/

struct user_coul_corr_data{
  double dA1;
  double dA2;
};

/*##############################################################################
// Prototype for user-supplied function.
//############################################################################*/

double
my_coulomb_correction(
  Libnucnet__Species *,
  double,
  double,
  double,
  struct user_coul_corr_data *
);

#ifdef __cplusplus
} /* extern "C" */
#endif
