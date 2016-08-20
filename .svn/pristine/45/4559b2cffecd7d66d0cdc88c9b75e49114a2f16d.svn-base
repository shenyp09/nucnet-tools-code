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
//! \brief Code for user-defined screening functions.
//!
////////////////////////////////////////////////////////////////////////////////

#include "screen.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */
namespace user
{

//##############################################################################
// set_screening_function().
//##############################################################################

void
set_screening_function( nnt::Zone& zone )
{

  Libnucnet__Zone__setScreeningFunction(
    zone.getNucnetZone(),
    (Libnucnet__Zone__screeningFunction) screening_function,
    NULL
  );

  zone.updateFunction(
    nnt::s_SCREENING_DATA_FUNCTION,
    static_cast<boost::function<boost::any()> >(
      boost::bind( get_screening_data, boost::ref( zone ) )
    ),
    "The data retriever for the default screening function."
  );

}

//##############################################################################
// get_screening_data.
//##############################################################################

screening_data_t
get_screening_data( nnt::Zone& zone )
{

  screening_data_t screening_data( zone );

  return screening_data;

}

//##############################################################################
// screening_function().
//##############################################################################

void
screening_function(
  Libnucnet__Zone * p_zone,
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  double d_rho,
  double d_Ye,
  double * p_forward,
  double * p_reverse
)
{

  Libnucnet__Species__nseCorrectionFactorFunction pf = NULL;
  double d_screen_f, d_screen_r, d_nse_corr = 1;

  reaction_screening_function(
    Libnucnet__Zone__getNet( p_zone ),
    p_reaction,
    d_t9,
    d_rho,
    d_Ye,
    Libnucnet__Zone__getScreeningData( p_zone ),
    &d_screen_f,
    &d_screen_r
  );

  pf = Libnucnet__Zone__getNseCorrectionFactorFunction( p_zone );

  if( pf )
  {
    d_nse_corr =
      Libnucnet__Net__computeReverseRatioCorrectionFactorForReaction(
        Libnucnet__Zone__getNet( p_zone ),
        p_reaction,
        d_t9,
        d_rho,
        d_Ye,
        pf,
        Libnucnet__Zone__getNseCorrectionFactorData( p_zone )
      );
  }

  if( d_screen_f >= d_screen_r )
  {
    *p_forward *= d_screen_f;
    *p_reverse *= d_screen_f * d_nse_corr;
  }
  else
  {
    *p_forward *= d_screen_r / d_nse_corr;
    *p_reverse *= d_screen_r;
  }

}

//##############################################################################
// reaction_screening_function().
//##############################################################################

void
reaction_screening_function(
  Libnucnet__Net * p_net,
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  double d_rho,
  double d_Ye,
  void *p_user_data,
  double *p_screen_f,
  double *p_screen_r
)
{

  unsigned int i_z1, i_a1;

  //============================================================================
  // Compute the screening in the forward direction.
  //============================================================================

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list( p_reaction );

  i_z1 = 0;
  i_a1 = 0;
  *p_screen_f = 1;

  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
  {

    unsigned int i_z2 =
      Libnucnet__Species__getZ(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( p_net ),
          Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        )
      ); 
    
    unsigned int i_a2 =
      Libnucnet__Species__getA(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( p_net ),
          Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        )
      ); 
    
    if( i_z1 != 0 && i_a1 != 0 )
    {
      *p_screen_f *=
        pair_screening_function(
          d_t9,
          d_rho,
          d_Ye,
          i_z1,
          i_a1,
          i_z2,
          i_a2,
          p_user_data
        );
    }

    i_z1 += i_z2;
    i_a1 += i_a2;

  }

  //============================================================================
  // Compute the screening in the reverse direction.
  //============================================================================

  nnt::reaction_element_list_t product_list =
    nnt::make_reaction_nuclide_product_list( p_reaction );

  i_z1 = 0;
  i_a1 = 0;
  *p_screen_r = 1;

  BOOST_FOREACH( nnt::ReactionElement product, product_list )
  {

    unsigned int i_z2 =
      Libnucnet__Species__getZ(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( p_net ),
          Libnucnet__Reaction__Element__getName(
            product.getNucnetReactionElement()
          )
        )
      ); 
    
    unsigned int i_a2 =
      Libnucnet__Species__getA(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( p_net ),
          Libnucnet__Reaction__Element__getName(
            product.getNucnetReactionElement()
          )
        )
      ); 
    
    if( i_z1 != 0 && i_a1 != 0 )
    {
      *p_screen_r *=
        pair_screening_function(
          d_t9,
          d_rho,
          d_Ye,
          i_z1,
          i_a1,
          i_z2,
          i_a2,
          p_user_data
        );
    }

    i_z1 += i_z2;
    i_a1 += i_a2;

  }

}

//##############################################################################
// Base screening function.
//##############################################################################

double
pair_screening_function(
  double d_T9,
  double d_rho,
  double d_Ye,
  unsigned int i_Z1,
  unsigned int i_A1,
  unsigned int i_Z2,
  unsigned int i_A2,
  void *p_data
)
{

  screening_data_t screening_data =
    boost::any_cast<screening_data_t>( *(boost::any *) p_data );

  double d_A_bar;

  double d_Gamma_e;
  double d_Gamma_effective;

  double d_H_strong;
  double d_H_weak;
  double d_H_intermediate;

  double d_f = 0.;
  double d_Ye2 = screening_data.dYe2;

  d_A_bar = (double) i_A1 * (double) i_A2 / ( (double) i_A1 + (double) i_A2 ); 

  d_Gamma_e = calculate_gamma_e( d_T9,d_rho,d_Ye );
  d_Gamma_effective = calculate_gamma_effective( i_Z1,i_Z2, d_Gamma_e );

  if( d_Gamma_effective < 0.3 )
  {
    d_H_weak =
      weak_screening_factor( i_Z1, i_Z2, d_A_bar, d_T9, d_rho, d_Ye, d_Ye2);
    d_f = exp(d_H_weak);
  }

  if( (d_Gamma_effective >= 0.3) && (d_Gamma_effective < 0.8) )
  {
    d_H_weak =
      weak_screening_factor( i_Z1, i_Z2, d_A_bar, d_T9, d_rho, d_Ye, d_Ye2);
    d_H_strong =
      strong_screening_factor(
        (double) i_Z1,
        (double) i_A1,
        (double) i_Z2,
        (double) i_A2,
        d_T9,
        d_Gamma_e,
        d_Gamma_effective
      );
    d_H_intermediate = intermediate_screening_factor( d_H_weak, d_H_strong);
    d_f = exp(d_H_intermediate);
  }

  if( (d_Gamma_effective >= 0.8) && (d_Gamma_effective <= 168) )
  {
    d_H_strong =
      strong_screening_factor(
        (double) i_Z1,
        (double) i_A1,
        (double) i_Z2,
        (double) i_A2,
        d_T9,
        d_Gamma_e,
        d_Gamma_effective
      );
    d_f = exp(d_H_strong);
  }

  if( d_Gamma_effective > 168 ) {
    LIBNUCNET__ERROR(
      "This is beyond the strong screening regime and not covered by this approximation."
    );
  }

  return d_f;

}

//##############################################################################
// intermediate_screening_factor
//##############################################################################

double intermediate_screening_factor( double d_H_weak, double d_H_strong ) {

//==============================================================================
//
// H_{intermediate} =
//    (H_{weak} H_{strong}) / (H_{weak}^{2} + H_{strong}^2)^{1/2}
//
//==============================================================================

  double d_H_intermediate;
  double d_denominator;

  d_denominator = pow(d_H_weak,2.0) + pow(d_H_strong, 2.0);
  d_denominator = pow(d_denominator,1.0/2.0);

  d_H_intermediate =  d_H_weak*d_H_strong/d_denominator;

  return d_H_intermediate;

}


//##############################################################################
// weak_screening_factor
//##############################################################################

double
weak_screening_factor(
  unsigned int i_Z1,
  unsigned int i_Z2,
  double d_A_bar,
  double d_T9,
  double d_rho,
  double d_Ye,
  double d_Ye2
)
{

//=============================================================================
//
// H_{weak} = Z_{1} Z_{2} \bar{A}^{-1/2} (Y_{e,2} + \theta_{e}Y_{e})^{1/2}
//  (1.88x10^{8}) (\rho/(T^{3})^{1/2}
//
// Y_{e,2} = \sum_{i} Z_{i}^{2} Y_{i}
//
// **WARNINGS**
//
//  1) Equation (A14) in Wallace (1982) contains errors. The above has these
//  errors corrected.
//  2) The above includes an electron degeneracy factor which would make the
//     Y_{e} term go to zero in degenerate matter. 
// 
//
//==============================================================================

 double d_H_weak;

  d_H_weak = i_Z1*i_Z2*1.88E8;
  d_H_weak = d_H_weak * pow(d_A_bar, -1.0/2.0);
  d_H_weak = d_H_weak * pow(D_THETA_E*d_Ye + d_Ye2, 1.0/2.0);
  d_H_weak = d_H_weak * pow(d_rho, 1.0/2.0) * pow(d_T9*1.0E9, -3.0/2.0);

  return d_H_weak;

}


//##############################################################################
// strong_screening_factor
//############################################################################*/

double
strong_screening_factor(
  double Z1,
  double A1,
  double Z2,
  double A2,
  double d_T9,
  double d_Gamma_e,
  double d_Gamma_effective
)
{

//=============================================================================
//
// \tau_{1,2} = 4.24872((Z_{1}^{2} Z_{2}^{2} \hat{A}_{1,2})/(T_{9}))^{1/3}
//
// \hat{A}_{1,2} = A_{1}A_{2}\over{A_{1} + A_{2}}
//
//
// C = 0.896434\Gamma_{effective}\tilde{z} -
//     3.44740(\Gamma_{effective})^{1/4}\tilde{\varsigma} -
//     0.5551(\ln{\Gamma_{e}} + (5/3)\ln{[(Z_{1}Z_{2}/(Z_{1} + Z_{2})]}) - 2.996
//
// \tilde{z} = (Z_{1} + Z_{2})^{5/3} - Z_{1}^{5/3} - Z_{2}^{5/3}
//
// \tilde{\varsigma} = (Z_{1} + Z_{2})^{5/12} - Z_{1}^{5/12} - Z_{2}^{5/12}
//
//
// b = 3.0\Gamma_{effective} / \tau_{1,2}
//
//
// H_{strong} = C - \tau_{1,2}/3((5/32)b^{3} - 0.014b^{4} -0.128b^{5}) - 
//             \Gamma_{effective}(0.0055b^{4} - 0.0098b^{5} + 0.0048b^{6})
//
//==============================================================================

  double d_A_reduced;
  double d_tau12;

  double d_z53;
  double d_z512;

  double d_C;

  double d_b;

  double d_H_strong;

  d_A_reduced = ( A1* A2 )/( A1 + A2);

  d_tau12 = gsl_pow_2( Z1 * Z1 ) * d_A_reduced/d_T9;
  d_tau12 = pow(d_tau12, 1.0/3.0);
  d_tau12 = 4.24872 * d_tau12;

  d_z53 =
    pow( Z1 + Z2, 5.0/3.0 ) -
    pow( Z1, 5.0/3.0) - pow( Z2, 5.0/3.0);
  d_z512 = pow( Z1 + Z2, 5.0/12.0) -
    pow( Z1, 5.0/12.0) - pow( Z2, 5.0/12.0);

  d_C = 0.896434 * d_Gamma_e * d_z53;
  d_C -= 3.44740 * pow(d_Gamma_e,1.0/4.0) * d_z512;
  d_C -= 0.5551 * log(d_Gamma_e);
  d_C -= 0.5551 * (5.0/3.0) * log( Z1 * Z2)/ ( Z1 + Z2 );
  d_C -= 2.996;

  d_b = 3.0 * d_Gamma_effective / d_tau12;

  d_H_strong =
    d_C - (d_tau12/3.0)*((5.0/32.0)*pow(d_b, 3.0) - 0.014*pow(d_b,4.0) -
    0.128*pow(d_b,5.0));

  d_H_strong += - d_Gamma_effective*(0.0055*pow(d_b,4.0) -
     0.0098*pow(d_b,5.0) + 0.0048*pow(d_b,6.0));

  return d_H_strong; 

}

//##############################################################################
// calculate_gamma_e
//##############################################################################

double calculate_gamma_e ( double d_T9, double d_rho, double d_Ye)
{

//=============================================================================
//
// \Gamma_{e} = (e^2/a_e KT)
//
// a_e = ( 4 * \pi * \rho * N_A * Y_e / 3 )^{1/3}
//
// a_e is the electron cloud radius.  Notice that gsl gives the electron
// charge e in abamps-s.  To convert to esu, multiply by the speed of light.
//==============================================================================

  double d_a_e;   // The electron cloud radius.
  double d_Gamma_e;

  d_a_e = 4.0*M_PI*d_rho*GSL_CONST_NUM_AVOGADRO*d_Ye/3.0;
  d_a_e = pow(d_a_e, -1.0/3.0);

  d_Gamma_e =
    gsl_pow_2(
      GSL_CONST_CGSM_ELECTRON_CHARGE * GSL_CONST_CGSM_SPEED_OF_LIGHT
    ) /
    ( d_a_e * GSL_CONST_CGSM_BOLTZMANN * d_T9 * GSL_CONST_NUM_GIGA );

  return d_Gamma_e;

}

//##############################################################################
// calculate_gamma_effective
//##############################################################################

double
calculate_gamma_effective(
  unsigned int i_Z1,
  unsigned int i_Z2,
  double d_Gamma_e
)
{

//==============================================================================
//
// \Gamma_{effective} = (2 {\over Z_{1} + Z_{2}})^{1/3} Z_{1} Z_{2} \Gamma_{e}
//
//==============================================================================

  double d_Gamma_effective;

  d_Gamma_effective = 2.0 / (i_Z1 + i_Z2);
  d_Gamma_effective = pow(d_Gamma_effective,1.0/3.0);
  d_Gamma_effective = d_Gamma_effective * i_Z1 * i_Z2 * d_Gamma_e;

  return d_Gamma_effective;

}

} // namespace user
