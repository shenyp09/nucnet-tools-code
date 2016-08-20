////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2015-2016 Clemson University.
//
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
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file hydro_helper.h
//! \brief A header file to define useful hydro helper routines.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef HYDRO_HELPER_H
#define HYDRO_HELPER_H

#include <algorithm>

#include <Libnucnet.h>

#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/assign.hpp>
#include <algorithm>

#include "nnt/iter.h"
#include "nnt/string_defs.h"

#include "user/evolve.h"

/**
 * @brief A namespace for user-defined functions.
 */
namespace user
{

typedef std::vector< double > state_type;

//##############################################################################
// Prototypes.
//##############################################################################

double acceleration( nnt::Zone&, const state_type&, const double );

double rho_function( nnt::Zone&, const state_type& );

double d_ln_rho_dt_function( nnt::Zone&, const state_type& );

void evolve_function( nnt::Zone&, Libnucnet__NetView *, const double );
 
double t9_from_entropy_root( double, nnt::Zone&, Libnucnet__NetView * );

double compute_entropy( nnt::Zone& );

double t9_function( nnt::Zone&, Libnucnet__NetView * );

void
observer_function(
  nnt::Zone&,
  const state_type&,
  const state_type&,
  const double
);

} // namespace user

#endif // HYDRO_HELPER_H
