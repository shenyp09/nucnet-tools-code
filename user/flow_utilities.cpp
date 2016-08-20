////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
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
//! \file
//! \brief Code for routines related to reaction flows.
//!
////////////////////////////////////////////////////////////////////////////////

#include "user/flow_utilities.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */
namespace user
{

//############################################################################
// make_flow_data_tuple().
//##########################################################################//

/**
 * \brief Construct a flow data tuple from zone data.
 *
 * \param zone A Nucnet Tools zone.
 * \return A tuple containing the zone Ye, screening function data, screening
 *         applier data, and nse correction factor function data.
*/

flow_data_tuple_t
make_flow_data_tuple(
  nnt::Zone& zone
)
{

  boost::any screening_data = 0;
  boost::any nse_correction_data = 0;

  if( Libnucnet__Zone__getScreeningFunction( zone.getNucnetZone() ) )
  {
    screening_data =
      boost::any_cast<boost::function<boost::any()> >(
        zone.getFunction( nnt::s_SCREENING_DATA_FUNCTION )
      )();
  }

  if( Libnucnet__Zone__getNseCorrectionFactorFunction( zone.getNucnetZone() ) )
  {
    nse_correction_data =
      boost::any_cast<boost::function<boost::any()> >(
        zone.getFunction( nnt::s_NSE_CORRECTION_FACTOR_DATA_FUNCTION )
      )();
  }

  return
    boost::make_tuple(
      user::compute_cluster_abundance_moment( zone, "", "z", 1. ),
      screening_data,
      nse_correction_data
    );

}

//############################################################################
// compute_rates_for_reaction_in_zone().
//##########################################################################//

/**
 * \brief Compute the forward and reverse rates for a reaction in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_reaction A pointer to the Libnucnet__Reaction.
 * \param data_tuple A tuple containing data for the screening and reverse
 *                      ratio functions (optional).
 * \return A pair in which the first element is the forward rate and the
 *         second is the reverse.
 */

std::pair<double,double>
compute_rates_for_reaction_in_zone(
  nnt::Zone &zone,
  Libnucnet__Reaction *p_reaction,
  flow_data_tuple_t& data_tuple
)
{

  double d_forward_rate, d_reverse_rate;
  Libnucnet__Zone__screeningFunction pfFunc = NULL;

  //==========================================================================
  // Store data.
  //==========================================================================

  double d_T9 = zone.getProperty<double>( nnt::s_T9 );
  double d_rho = zone.getProperty<double>( nnt::s_RHO );
  double d_Ye = data_tuple.get<I_TUPLE_YE>();

  //==========================================================================
  // Compute rates.
  //==========================================================================

  Libnucnet__Net__computeRatesForReaction(
    Libnucnet__Zone__getNet( zone.getNucnetZone() ),
    p_reaction,
    d_T9,
    d_rho,
    Libnucnet__Zone__getDataForUserRateFunction(
      zone.getNucnetZone(),
      Libnucnet__Reaction__getRateFunctionKey( p_reaction )
    ),
    &d_forward_rate,
    &d_reverse_rate
  );

  if( Libnucnet__Reaction__isWeak( p_reaction ) )
  {

    d_reverse_rate =
      compute_reverse_weak_rate_for_reaction(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
        p_reaction,
        d_forward_rate,
        d_T9,
        d_rho,
        zone.getProperty<double>( nnt::s_MUEKT ),
        compute_thermo_quantity(
          zone,
          nnt::s_CHEMICAL_POTENTIAL_KT,
          nnt::s_NEUTRINO_E
        )
      );

  }

  pfFunc = Libnucnet__Zone__getScreeningFunction( zone.getNucnetZone() );

  if( pfFunc )
  {

    Libnucnet__Zone__setScreeningFunction(
      zone.getNucnetZone(),
      (Libnucnet__Zone__screeningFunction) pfFunc,
      &(data_tuple.get<I_TUPLE_SCREENING_DATA>())
    );

    Libnucnet__Zone__setNseCorrectionFactorFunction(
      zone.getNucnetZone(),
      (Libnucnet__Species__nseCorrectionFactorFunction)
        Libnucnet__Zone__getNseCorrectionFactorFunction(
          zone.getNucnetZone()
        ),
      &(data_tuple.get<I_TUPLE_NSE_CORR_DATA>())
    );
  
    pfFunc(
      zone.getNucnetZone(),
      p_reaction,
      d_T9,
      d_rho,
      d_Ye,
      &d_forward_rate,
      &d_reverse_rate
    );

  }
  
  //==========================================================================
  // Modify rates for reaction.
  //==========================================================================

  modify_rates_for_reaction(
    zone,
    p_reaction,
    d_forward_rate,
    d_reverse_rate
  );

  //==========================================================================
  // Return.
  //==========================================================================

  return std::make_pair( d_forward_rate, d_reverse_rate );

}

std::pair<double,double>
compute_rates_for_reaction_in_zone(
  nnt::Zone &zone,
  Libnucnet__Reaction *p_reaction
)
{

  flow_data_tuple_t data_tuple = make_flow_data_tuple( zone );

  return
    compute_rates_for_reaction_in_zone(
      zone,
      p_reaction,
      data_tuple
    );

}

//############################################################################
// compute_reaction_abundance_product_in_zone().
//##########################################################################//

/**
 * \brief Compute the product of abundances for a reaction in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param element_list A NucNet Tools reaction element list.
 * \param p_species_to_exclude A species whose abundance should be excluded
 *        once from the product (optional).
 * \return A double giving the product.
 */

double
compute_reaction_abundance_product_in_zone(
  nnt::Zone& zone,
  nnt::reaction_element_list_t& element_list,
  Libnucnet__Species * p_species_to_exclude
)
{

  double d_result = 1;
  Libnucnet__Species * p_exclude = p_species_to_exclude;

  BOOST_FOREACH( nnt::ReactionElement element, element_list )
  {

    Libnucnet__Species * p_species =
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
        Libnucnet__Reaction__Element__getName(
          element.getNucnetReactionElement()
        )
      );

    if( p_species != p_exclude )
    {
      d_result *=
        Libnucnet__Zone__getSpeciesAbundance(
    	  zone.getNucnetZone(),
          p_species
        );
    }
    else
    {
      p_exclude = NULL;   // To exclude species only once.
    }

  }

  return d_result;

}

double
compute_reaction_abundance_product_in_zone(
  nnt::Zone& zone,
  nnt::reaction_element_list_t& element_list
)
{

  Libnucnet__Species * p_null = NULL;

  return
    compute_reaction_abundance_product_in_zone( zone, element_list, p_null );

}

//############################################################################
// compute_flows_for_reaction().
//##########################################################################//

/**
 * \brief Compute the forward and reverse flows for a reaction in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_reaction A pointer to the Libnucnet__Reaction.
 * \param data_tuple A tuple containing data for the screening and reverse
 *                      ratio functions (optional).
 * \return A pair in which the first element is the forward flow and the
 *         second is the reverse.
 */

std::pair<double,double>
compute_flows_for_reaction(
  nnt::Zone &zone,
  Libnucnet__Reaction *p_reaction,
  flow_data_tuple_t& data_tuple
)
{

  double d_forward_rate, d_reverse_rate;
  double d_f, d_r;
  size_t i_elements;

  //==========================================================================
  // Compute rates.
  //==========================================================================

  boost::tie( d_forward_rate, d_reverse_rate ) =
    compute_rates_for_reaction_in_zone( zone, p_reaction, data_tuple );

  //==========================================================================
  // Forward flow.
  //==========================================================================

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list( p_reaction );

  i_elements = reactant_list.size();

  d_f =
    d_forward_rate
    *
    pow( zone.getProperty<double>( nnt::s_RHO ), (double) i_elements - 1. )
    *
    compute_reaction_abundance_product_in_zone( zone, reactant_list ) 
    /
    Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

  //==========================================================================
  // Reverse flow.
  //==========================================================================

  nnt::reaction_element_list_t product_list =
    nnt::make_reaction_nuclide_product_list( p_reaction );

  i_elements = product_list.size();

  d_r =
    d_reverse_rate
    *
    pow( zone.getProperty<double>( nnt::s_RHO ), (double) i_elements - 1. )
    *
    compute_reaction_abundance_product_in_zone( zone, product_list )
    /
    Libnucnet__Reaction__getDuplicateProductFactor( p_reaction );

  //==========================================================================
  // Return the flows.
  //==========================================================================
  
  return std::make_pair( d_f, d_r ); 

}

std::pair<double,double>
compute_flows_for_reaction(
  nnt::Zone &zone,
  Libnucnet__Reaction *p_reaction
)
{

  flow_data_tuple_t data_tuple = make_flow_data_tuple( zone );

  return
    compute_flows_for_reaction(
      zone,
      p_reaction,
      data_tuple
    );

}

//############################################################################
// compute_forward_flow_vector().
//##########################################################################//

/**
 * \brief Compute the forward flow vector for a set of nuclei and reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions to use.
 *         the network nuclei.
 * \param s_reac_xpath A string giving an XPath expression restricting
 *         the network reactions.
 * \return A pointer to a new gsl vector containing the flows for each species.
 */

gsl_vector *
compute_forward_flow_vector(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  gsl_vector *p_forward_flow;
  std::pair<double,double> flows;

  p_forward_flow =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );

  flow_data_tuple_t flow_data_tuple = make_flow_data_tuple( zone );
    
  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction(),
        flow_data_tuple
      );

    nnt::reaction_element_list_t element_list =
      nnt::make_reaction_nuclide_reactant_list( reaction.getNucnetReaction() );

    BOOST_FOREACH( nnt::ReactionElement element, element_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_forward_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_forward_flow,
          Libnucnet__Species__getIndex( p_species )
        ) + flows.first
      );

    }

  }

  return p_forward_flow;

}

//############################################################################
// compute_reverse_flow_vector().
//##########################################################################//

/**
 * \brief Compute the reverse flow vector for a set of nuclei and reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions to use.
 * \return A pointer to a new gsl vector containing the flows for each species.
 */

gsl_vector *
compute_reverse_flow_vector(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  gsl_vector * p_reverse_flow;
  std::pair<double,double> flows;

  p_reverse_flow =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    );

  flow_data_tuple_t flow_data_tuple = make_flow_data_tuple( zone );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );
    
  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction(),
        flow_data_tuple
      );

    nnt::reaction_element_list_t element_list =
      nnt::make_reaction_nuclide_product_list( reaction.getNucnetReaction() );

    BOOST_FOREACH( nnt::ReactionElement element, element_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_reverse_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_reverse_flow,
          Libnucnet__Species__getIndex( p_species )
        ) + flows.second
      );

    }

  }

  return p_reverse_flow;

}

//############################################################################
// compute_flow_vector().
//##########################################################################//

/**
 * \brief Compute the net flow vector for a set of nuclei and reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions to use.
 */

gsl_vector *
compute_flow_vector(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  gsl_vector *p_flow;
  std::pair<double,double> flows;

  p_flow =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    );

  flow_data_tuple_t flow_data_tuple = make_flow_data_tuple( zone );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );
    
  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction(),
        flow_data_tuple
      );

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list(
        reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_flow,
          Libnucnet__Species__getIndex( p_species )
        ) + flows.first - flows.second
      );

    }

    nnt::reaction_element_list_t product_list =
      nnt::make_reaction_nuclide_product_list(
        reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement product, product_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            product.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_flow,
          Libnucnet__Species__getIndex( p_species )
        ) - flows.first + flows.second
      );

    }

  }

  return p_flow;

}

//############################################################################
// compute_total_flow().
//##########################################################################//

/**
 * \brief Compute the total net flow for a set of nuclei and their reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param s_type_flow A string giving the type of flow ("forward", "reverse",
 *    or "net").
 * \param A network view specifying the nuclei and reactions to use.
 * \return A double giving the sum of the net flows for the chosen reactions.
 */

double
compute_total_flow(
  nnt::Zone& zone,
  const char * s_type_flow,
  Libnucnet__NetView * p_net_view
)
{

  gsl_vector * p_flow;
  double d_result = 0;
  size_t i;

  if( strcmp( s_type_flow, nnt::s_FORWARD_FLOW ) == 0 )
    p_flow = compute_forward_flow_vector( zone, p_net_view );
  else if( strcmp( s_type_flow, nnt::s_REVERSE_FLOW ) == 0 )
    p_flow = compute_reverse_flow_vector( zone, p_net_view );
  else if( strcmp( s_type_flow, nnt::s_NET_FLOW ) == 0 )
    p_flow = compute_flow_vector( zone, p_net_view );
  else
  {
    std::cerr << "No such flow type in compute_total_flow." << std::endl;
    exit( EXIT_FAILURE );
  }

  for( i = 0; i < p_flow->size; i++ )
    d_result += gsl_vector_get( p_flow, i );

  gsl_vector_free( p_flow );

  return d_result;

}

//##############################################################################
// compute_zone_reactions_entropy_generation_rate_per_nucleon().
//##############################################################################

/**
 * \brief Compute the entropy generation rate per nucleon for reactions in
 *        a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A vector of pairs giving the reaction string and rate.
 */

std::vector<std::pair<std::string, double> >
compute_zone_reactions_entropy_generation_rate_per_nucleon(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  std::pair<double,double> flows;
  int i_delta_z = 0;
  double d_reaction_entropy_change, d_abund;
  double d_mu_nue_kT;
  std::vector<double> nse_corr_vector;
  size_t i = 0;
  Libnucnet__Species__nseCorrectionFactorFunction pfFunc = NULL;
 
  std::vector<std::pair<std::string,double> > result(
    Libnucnet__Reac__getNumberOfReactions(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    )
  );

  zone.updateProperty(
    nnt::s_MUEKT,
    compute_thermo_quantity(
      zone,
      nnt::s_CHEMICAL_POTENTIAL_KT,
      nnt::s_ELECTRON
    )
  );

  flow_data_tuple_t flow_data_tuple = make_flow_data_tuple( zone );

  pfFunc =
    Libnucnet__Zone__getNseCorrectionFactorFunction(
      zone.getNucnetZone()
    );

  if( pfFunc )
  {

    nnt::species_list_t species_list =
      nnt::make_species_list(
        Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) )
      );


    BOOST_FOREACH( nnt::Species species, species_list )
    {

      nse_corr_vector.push_back(
        pfFunc(
          species.getNucnetSpecies(),
          zone.getProperty<double>( nnt::s_T9 ),
          zone.getProperty<double>( nnt::s_RHO ),
          flow_data_tuple.get<I_TUPLE_YE>(),
          &(flow_data_tuple.get<I_TUPLE_NSE_CORR_DATA>())
        )
      );

    }

  }

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_net_view )
      )
    );

  d_mu_nue_kT =
    compute_thermo_quantity(
      zone,
      nnt::s_CHEMICAL_POTENTIAL_KT,
      nnt::s_NEUTRINO_E
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction(),
        flow_data_tuple
      );

    i_delta_z = 0;

    if( flows.first > 1.e-100 )
    {

      d_reaction_entropy_change =
	nnt::compute_reaction_nuclear_Qvalue(
	  Libnucnet__NetView__getNet( p_net_view ),
	  reaction.getNucnetReaction(),
	  nnt::d_ELECTRON_MASS_IN_MEV
	) /
        nnt::compute_kT_in_MeV( zone.getProperty<double>( nnt::s_T9 ) );

      nnt::reaction_element_list_t reactant_list =
	nnt::make_reaction_nuclide_reactant_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
      {

	p_species =
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) ),
	    Libnucnet__Reaction__Element__getName(
	      reactant.getNucnetReactionElement()
	    )
	  );

	i_delta_z += (int) Libnucnet__Species__getZ( p_species );

	d_abund =
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    p_species
	  );

	if( d_abund > 1.e-100 )
	{
	  d_reaction_entropy_change +=
	    log(
	      d_abund
	      /
	      Libnucnet__Species__computeQuantumAbundance(
		p_species,
		zone.getProperty<double>( nnt::s_T9 ),
		zone.getProperty<double>( nnt::s_RHO )
	      )
	    );
	}

        if( pfFunc )
        {
          d_reaction_entropy_change -=
            nse_corr_vector[ Libnucnet__Species__getIndex( p_species )];
        }

      }

      nnt::reaction_element_list_t product_list =
	nnt::make_reaction_nuclide_product_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement product, product_list )
      {

	p_species =
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) ),
	    Libnucnet__Reaction__Element__getName(
	      product.getNucnetReactionElement()
	    )
	  );

	i_delta_z -= (int) Libnucnet__Species__getZ( p_species );

	d_abund =
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    p_species
	  );

	if( d_abund > 1.e-100 )
	{
	  d_reaction_entropy_change -=
	    log(
	      d_abund
	      /
	      Libnucnet__Species__computeQuantumAbundance(
		p_species,
		zone.getProperty<double>( nnt::s_T9 ),
		zone.getProperty<double>( nnt::s_RHO )
	      )
	    );
	}

        if( pfFunc )
        {
          d_reaction_entropy_change +=
            nse_corr_vector[ Libnucnet__Species__getIndex( p_species )];
        }

      }

      if( i_delta_z != 0 )
      {

	d_reaction_entropy_change +=
	  (double) i_delta_z *
	  (
	     zone.getProperty<double>( nnt::s_MUEKT ) +
	     nnt::d_ELECTRON_MASS_IN_MEV /
	     nnt::compute_kT_in_MeV( zone.getProperty<double>( nnt::s_T9 ) )
	  );

	if( d_mu_nue_kT != GSL_NEGINF )
	{
	  d_reaction_entropy_change -= (double) i_delta_z * d_mu_nue_kT;
	}

      }

      result[i++] =
        std::make_pair(
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
          d_reaction_entropy_change * ( flows.first - flows.second )
        );

    }

  }

  return result;

}
           
//##############################################################################
// compute_entropy_generation_rate().
//##############################################################################

/**
 * \brief Compute the entropy generation rate per nucleon in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A double giving the entropy change rate per nucleon.
 */

double
compute_entropy_generation_rate(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  double d_result = 0;

  std::vector<std::pair<std::string, double> > result =
    compute_zone_reactions_entropy_generation_rate_per_nucleon(
      zone,
      p_net_view
    );

  for( size_t i = 0; i < result.size(); i++ )
  {
    d_result += result[i].second;
  }

  return d_result;

}

//##############################################################################
// compute_energy_generation_rate_per_nucleon_for_reaction().
//##############################################################################

double
compute_energy_generation_rate_per_nucleon_for_reaction(
  nnt::Zone& zone,
  Libnucnet__Reaction * p_reaction,
  double d_electron_energy_part,
  flow_data_tuple_t flow_data_tuple
)
{

  int i_y = 0, i_z = 0;
  std::pair<double,double> flows;

  flows =
    compute_flows_for_reaction( zone, p_reaction, flow_data_tuple );

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list(
      p_reaction
    );
        
  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
  {

    i_z -=
      (int)
      Libnucnet__Species__getZ(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        )
      );

    ++i_y;

  }

  nnt::reaction_element_list_t product_list =
    nnt::make_reaction_nuclide_product_list( p_reaction );
        
  BOOST_FOREACH( nnt::ReactionElement product, product_list )
  {

    i_z +=
      (int)
      Libnucnet__Species__getZ(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            product.getNucnetReactionElement()
          )
        )
      );

    --i_y;

  }

  return
    (
      (
        nnt::compute_reaction_nuclear_Qvalue(
          Libnucnet__Zone__getNet( zone.getNucnetZone() ),
          p_reaction,
          nnt::d_ELECTRON_MASS_IN_MEV
        ) *
        GSL_CONST_CGSM_ELECTRON_VOLT *
        GSL_CONST_NUM_MEGA
      )
      +
      i_y *
      (
        1.5 *
        GSL_CONST_CGSM_BOLTZMANN *
        zone.getProperty<double>( nnt::s_T9 ) *
        GSL_CONST_NUM_GIGA
      )
      +
      i_z * 
      (
        d_electron_energy_part
        +
        nnt::d_ELECTRON_MASS_IN_MEV *
        GSL_CONST_CGSM_ELECTRON_VOLT *
        GSL_CONST_NUM_MEGA
      )
    ) * ( flows.first - flows.second );
/*
        gsl_pow_2(
          zone.getProperty<double>( nnt::s_T9 ) *
          GSL_CONST_NUM_GIGA
        ) *
          GSL_CONST_CGSM_BOLTZMANN *
          compute_dlnG_dT(
            reaction.getNucnetSpecies(),
            zone.getProperty<double>( nnt::s_T9 ) *
              GSL_CONST_NUM_GIGA
          )
*/

}

//##############################################################################
// compute_energy_generation_rate_per_nucleon().
//##############################################################################

/**
 * \brief Compute the energy generation rate per nucleon in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A double giving the energy generation rate (ergs per second per
 *         nucleon).
 */

double
compute_energy_generation_rate_per_nucleon(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  double d_result = 0, d_electron_energy_part;

  zone.updateProperty(
    nnt::s_MUEKT,
    compute_thermo_quantity(
      zone,
      nnt::s_CHEMICAL_POTENTIAL_KT,
      nnt::s_ELECTRON
    )
  );

  d_electron_energy_part =
    GSL_CONST_CGSM_BOLTZMANN *
    gsl_pow_2(
      zone.getProperty<double>( nnt::s_T9 ) *
      GSL_CONST_NUM_GIGA
    ) *
    compute_thermo_quantity(
      zone,
      nnt::s_T_DERIVATIVE_CHEMICAL_POTENTIAL_KT,
      nnt::s_ELECTRON
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_net_view )
      )
    );

  flow_data_tuple_t flow_data_tuple = make_flow_data_tuple( zone );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {
    d_result +=
      compute_energy_generation_rate_per_nucleon_for_reaction(
        zone,
        reaction.getNucnetReaction(),
        d_electron_energy_part,
        flow_data_tuple
      );
  }

  return d_result;

}
           
//##############################################################################
// compute_energy_generation_rate_per_gram().
//##############################################################################

/**
 * \brief Compute the energy generation rate per gram in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A double giving the energy generation rate (ergs per second
 *         per gram).
 */

double
compute_energy_generation_rate_per_gram(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  return
    compute_energy_generation_rate_per_nucleon( zone, p_net_view ) *
    GSL_CONST_NUM_AVOGADRO;

}
   
//##############################################################################
// compute_zone_reactions_energy_generation_rate_per_nucleon().
//##############################################################################

/**
 * \brief Compute the energy generation rate (ergs per second) per nucleon for
 *        reactions in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A vector of pairs giving the reaction string and rate.
 */

std::vector<std::pair<std::string,double> >
compute_zone_reactions_energy_generation_rate_per_nucleon(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  std::vector<std::pair<std::string,double> > result(
    Libnucnet__Reac__getNumberOfReactions(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    )
  );

  double d_electron_energy_part =
    GSL_CONST_CGSM_BOLTZMANN *
    gsl_pow_2(
      zone.getProperty<double>( nnt::s_T9 ) *
      GSL_CONST_NUM_GIGA
    ) *
    compute_thermo_quantity(
      zone,
      nnt::s_T_DERIVATIVE_CHEMICAL_POTENTIAL_KT,
      nnt::s_ELECTRON
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_net_view )
      )
    );

  flow_data_tuple_t flow_data_tuple = make_flow_data_tuple( zone );

  size_t i = 0;
  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {
    result[i++] =
      std::make_pair(
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        compute_energy_generation_rate_per_nucleon_for_reaction(
          zone,
          reaction.getNucnetReaction(),
          d_electron_energy_part,
          flow_data_tuple
        )
      );
  }

  return result;

}
           
//##############################################################################
// compute_zone_reactions_energy_generation_rate_per_gram().
//##############################################################################

/**
 * \brief Compute the energy generation rate (ergs per second) per gram for
 *        reactions in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A vector of pairs giving the reaction string and rate.
 */

std::vector<std::pair<std::string, double> >
compute_zone_reactions_energy_generation_rate_per_gram(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  std::vector<std::pair<std::string, double> > result =
    compute_zone_reactions_energy_generation_rate_per_nucleon(
      zone,
      p_net_view
    );

  for( size_t i = 0; i < result.size(); i++ )
  {
    result[i].second *= GSL_CONST_NUM_AVOGADRO;
  }

  return result;

}
           
//############################################################################
// update_flow_currents().
//##########################################################################//

void
update_flow_currents( nnt::Zone& zone, nnt::Zone& flow_current_zone )
{

  if(
    zone.hasFunction( nnt::s_RATE_DATA_UPDATE_FUNCTION )
  )
  {

    boost::any_cast<boost::function<void( )> >(
      zone.getFunction( nnt::s_RATE_DATA_UPDATE_FUNCTION )
    )( );

  }

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet(
          zone.getNetView( "", "" )
        )
      )
    );

  double d_dt = zone.getProperty<double>( nnt::s_DTIME );

  flow_data_tuple_t flow_data_tuple = make_flow_data_tuple( zone );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    std::pair<double,double> flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction(),
        flow_data_tuple
      );

    double d_current = ( flows.first - flows.second ) * d_dt;

    if(
      flow_current_zone.hasProperty(
        nnt::s_FLOW_CURRENT, 
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
      )
    )
    {
      flow_current_zone.updateProperty(
        nnt::s_FLOW_CURRENT, 
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        flow_current_zone.getProperty<double>(
          nnt::s_FLOW_CURRENT,
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
        ) + d_current
      );
    }
    else
    {
      flow_current_zone.updateProperty(
        nnt::s_FLOW_CURRENT,
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        d_current
      );
    }

  }

}

} // namespace user
