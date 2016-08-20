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
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for routines related to reaction flows.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_FLOW_UTILITIES_H
#define USER_FLOW_UTILITIES_H

#include <iostream>
#include <utility>

#include <Libnucnet.h>
#include <gsl/gsl_vector.h>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

#include "thermo.h"
#include "user/rate_modifiers.h"
#include "screen.h"
#include "nse_corr.h"
#include "weak_utilities.h"

namespace user
{

enum
data_tuple_indices{
  I_TUPLE_YE,
  I_TUPLE_SCREENING_DATA,
  I_TUPLE_NSE_CORR_DATA
};

/**
 * A type consisting of a tuple containing the zone Ye, screening function
 * data, and nse correction factor function data.
 */
typedef
boost::tuple<double,boost::any,boost::any>
  flow_data_tuple_t;

flow_data_tuple_t
make_flow_data_tuple(
  nnt::Zone&
);

double
compute_reaction_abundance_product_in_zone(
  nnt::Zone &,
  nnt::reaction_element_list_t&,
  Libnucnet__Reaction *
);

double
compute_reaction_abundance_product_in_zone(
  nnt::Zone &,
  nnt::reaction_element_list_t&
);

std::pair<double,double>
compute_rates_for_reaction_in_zone(
  nnt::Zone &,
  Libnucnet__Reaction *,
  flow_data_tuple_t&
);

std::pair<double,double>
compute_rates_for_reaction_in_zone(
  nnt::Zone &,
  Libnucnet__Reaction *
);

std::pair<double,double>
compute_flows_for_reaction(
  nnt::Zone&,
  Libnucnet__Reaction *,
  flow_data_tuple_t&
);

std::pair<double,double>
compute_flows_for_reaction( nnt::Zone&, Libnucnet__Reaction * );

std::pair<double,double>
compute_modified_flows_for_reaction(
  nnt::Zone&,
  Libnucnet__Reaction *
);

gsl_vector *
compute_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_modified_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_forward_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_reverse_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_modfied_flow_vector( nnt::Zone&, Libnucnet__NetView * );

double
compute_total_flow(
  nnt::Zone&,
  const char *,
  Libnucnet__NetView *
);

std::vector<std::pair<std::string, double> >
compute_zone_reactions_entropy_generation_rate_per_nucleon(
  nnt::Zone&,
  Libnucnet__NetView *
);

double
compute_entropy_generation_rate(
  nnt::Zone&,
  Libnucnet__NetView *
);

double
compute_energy_generation_rate_per_nucleon(
  nnt::Zone&,
  Libnucnet__NetView *
);

double
compute_energy_generation_rate_per_gram(
  nnt::Zone&,
  Libnucnet__NetView *
);

std::vector<std::pair<std::string,double> >
compute_zone_reactions_energy_generation_rate_per_nucleon(
  nnt::Zone&,
  Libnucnet__NetView *
);

std::vector<std::pair<std::string,double> >
compute_zone_reactions_energy_generation_rate_per_gram(
  nnt::Zone&,
  Libnucnet__NetView *
);

void
update_flow_currents( nnt::Zone&, nnt::Zone& );

} // namespace user

#endif // USER_FLOW_UTILITIES_H

