//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer and George C.
// Jordan, IV.
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
//! \brief A header file defining screening functions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_SCREEN_H
#define USER_SCREEN_H

#include <Libnucnet.h>

#include "nnt/wrappers.hpp"
#include "nnt/string_defs.h"

#include "user/network_utilities.h"

namespace user
{

//##############################################################################
// Defines.
//##############################################################################

#define D_THETA_E   0.0   /* setting to 1.0 for non-degenerate matter.
                             Set to 0.0 for the completely degenerate case.
                          */

//##############################################################################
// User-supplied data structure.
//##############################################################################

typedef struct screening_data_t
{
  double dYe2;

  screening_data_t( nnt::Zone& zone ) 
  {
    dYe2 = user::compute_cluster_abundance_moment( zone, "", "z", 2. );
  }

} screening_data_t;

//##############################################################################
// Prototypes for user-supplied functions.
//##############################################################################

void
set_screening_function( nnt::Zone& );

screening_data_t
get_screening_data( nnt::Zone& );

void
screening_function(
  Libnucnet__Zone *,
  Libnucnet__Reaction *,
  double,
  double,
  double,
  double *,
  double *
);

void
reaction_screening_function(
  Libnucnet__Net *,
  Libnucnet__Reaction *,
  double,
  double,
  double,
  void *,
  double *,
  double *
);

double
pair_screening_function(
  double,
  double,
  double,
  unsigned int,
  unsigned int,
  unsigned int,
  unsigned int,
  void *
);

double calculate_gamma_e( double, double, double );

double calculate_gamma_effective( unsigned int, unsigned int, double );

double
strong_screening_factor(
  double,
  double,
  double,
  double,
  double,
  double,
  double
);

double
weak_screening_factor(
  unsigned int, unsigned int, double, double, double, double, double
);

double intermediate_screening_factor( double, double );

} // namespace user

#endif // USER_SCREEN_H
