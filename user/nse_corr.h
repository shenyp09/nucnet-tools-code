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
//! \brief A header for user-defined Coulomb correction functions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_NSE_CORR_H
#define USER_NSE_CORR_H

#include <iostream>

#include "nnt/wrappers.hpp"
#include "nnt/string_defs.h"

#include <Libnucnet.h>

#include "user/network_utilities.h"

namespace user
{

//##############################################################################
// User-supplied data structure.
//##############################################################################

typedef struct nse_corr_data_t
{

  double a, b, c, d, f1kT, o, beta, gamma;

  nse_corr_data_t( )
  {

    a = -0.898004;
    b =  0.9678;
    c =  0.22070;
    d = -0.86097;

    f1kT = -0.420;    // Slattery, Doolen, and DeWitt (Phys. Rev. A, 1980)
    o = a + 4. * (b - c) - f1kT;

    beta = a + b + c + d + sqrt( 3. ) / 2.;
    gamma = beta / ( ( 1. / sqrt( 3. ) ) + f1kT );

  }

} nse_corr_data_t;

//##############################################################################
// Prototypes for user-supplied functions.
//##############################################################################

void
set_nse_correction_function( nnt::Zone& );

nse_corr_data_t
get_nse_correction_data( );

double
nse_correction(
  Libnucnet__Species *, double, double, double, void *
);

double
species_coulomb_chemical_potential(
  Libnucnet__Species *, double, double, double, void *
);

double
species_coulomb_energy(
  Libnucnet__Species *, double, double, double, void *
);

double
species_coulomb_entropy(
  Libnucnet__Species *, double, double, double, void *
);

} // namespace user

#endif // USER_NSE_CORR_H
