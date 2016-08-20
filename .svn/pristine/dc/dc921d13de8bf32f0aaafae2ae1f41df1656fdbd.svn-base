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
//!
//! \file
//! \brief A header file to define routines and structures for user-defined
//!        rate functions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_RATE_FUNCTIONS_H
#define USER_RATE_FUNCTIONS_H

#include <iostream>
#include <string>

#include <Libnucnet.h>

#include "user/aa522a25.h"
#include "nnt/two_d_weak_rates.h"
#include "nnt/string_defs.h"

#include "user/neutrino_rate_functions.h"
#include "user/thermo.h"
#include "user/weak_utilities.h"

namespace user
{

//##############################################################################
// Prototypes.
//##############################################################################

void
set_rate_data_update_function( nnt::Zone& );

void
register_rate_functions( Libnucnet__Reac * );

void
update_two_d_weak_rate_functions_data(
  Libnucnet__Zone *,
  double,
  double,
  double,
  double
);

void
update_approximate_weak_rate_functions_data(
  Libnucnet__Zone *,
  double,
  double,
  double
);

void
update_rate_functions_data(
  nnt::Zone &
);

void
update_user_rate_functions_low_T_rates(
  Libnucnet__Reac *,
  const char *,
  const char *,
  double,
  double
);

} // namespace user

#endif // USER_RATE_FUNCTIONS_H
