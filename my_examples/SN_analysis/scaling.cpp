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
////////////////////////////////////////////////////////////////////////////////
 
////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Code for scaling flow arrows.
////////////////////////////////////////////////////////////////////////////////

#include "scaling.h"

//##############################################################################
// scale_graph_weights().
//##############################################################################

void
scale_graph_weights(
  my_graph_t& g,
  std::string s_type
)
{

  my_graph_t::edge_iterator ei, ei_end, e_next;

  //============================================================================
  // d_F is the factor in the logarithmic scaling for a flow arrow such
  // that the arrow width scales proportionally to
  //
  // 1 + (1 / d_F) * log10( Flow / Flow_max )
  //
  // Thus, for example, if a flow is 0.01 of the maximum flow, and if d_F = 5,
  // the arrow representing that flow will be 1 - 2 / 5 = 3 / 5.  Flows less
  // than 10^(-d_F) of the maximum are zeroed out.
  //============================================================================
  
  double d_F = 5;

  //============================================================================
  // Main code.
  //============================================================================
  
  double d_max = get_max_flow( g );

  if( GSL_SIGN( d_max ) == GSL_SIGN( -d_max ) ) d_max = 1.;

  boost::tie( ei, ei_end ) = boost::edges( g );

  //============================================================================
  // Choose scaling according to s_type.
  //============================================================================
  
  //----------------------------------------------------------------------------
  // Linear scaling.
  //----------------------------------------------------------------------------
  
  if( s_type == "linear" )
  {
    for( e_next = ei; ei != ei_end; ei = e_next )
    {
      g[*ei].setWeight( g[*ei].getWeight() / d_max );
      ++e_next;
      if( g[*ei].getWeight() < pow(10.,-d_F) ) boost::remove_edge( *ei, g );
    }
  }

  //----------------------------------------------------------------------------
  // Logarithmic scaling.
  //----------------------------------------------------------------------------
  
  else if( s_type == "logarithmic" )
  {
    for( e_next = ei; ei != ei_end; ei = e_next )
    {
      g[*ei].setWeight(
        1. + (1./d_F) * log10( g[*ei].getWeight() / d_max )
      );
      ++e_next;
      if( g[*ei].getWeight() < 0. ) boost::remove_edge( *ei, g );
    }
  }

  //----------------------------------------------------------------------------
  // Scaling not found.
  //----------------------------------------------------------------------------
  
  else
  {
    std::cerr << "No such scaling found." << std::endl;
    exit( EXIT_FAILURE );
  }

}

