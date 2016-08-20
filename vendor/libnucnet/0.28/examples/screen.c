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
//       Example of a user-supplied screening function.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "screen.h"

/*##############################################################################
// my_screening_function().
//############################################################################*/

void
my_screening_function(
  Libnucnet__Zone * p_zone,
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  double d_rho,
  double d_Ye,
  double * p_forward,
  double * p_reverse
)
{

  double d_screen_f, d_screen_r, d_nse_corr = 1;
  Libnucnet__Species__nseCorrectionFactorFunction pf = NULL;

  my_reaction_screening_function(
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

/*##############################################################################
// my_reaction_screening_function().
//############################################################################*/

void
my_reaction_screening_function(
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

  struct work_data wd;

  wd.pNet = p_net;
  wd.dT9 = d_t9;
  wd.dRho = d_rho;
  wd.dYe = d_Ye;
  wd.pUserData = p_user_data;

  wd.dScreen = 1;
  wd.iZ1 = 0;
  wd.iA1 = 0;

  Libnucnet__Reaction__iterateNuclideReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction)    
       compute_screening_factor_for_reaction,
    &wd
  );

  *p_screen_f = wd.dScreen;

  wd.dScreen = 1;
  wd.iZ1 = 0;
  wd.iA1 = 0;

  Libnucnet__Reaction__iterateNuclideProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction)    
       compute_screening_factor_for_reaction,
    &wd
  );

  *p_screen_r = wd.dScreen;

}

/*##############################################################################
// compute_screening_factor_for_reaction().
//############################################################################*/

int
compute_screening_factor_for_reaction(
  Libnucnet__Reaction__Element *p_element, void *p_data
)
{

  unsigned int i_z2, i_a2;

  struct work_data * p_wd = (struct work_data *) p_data;

  i_z2 =
    Libnucnet__Species__getZ(
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( p_wd->pNet ),
        Libnucnet__Reaction__Element__getName( p_element )
      )
    ); 
    
  i_a2 =
    Libnucnet__Species__getA(
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( p_wd->pNet ),
        Libnucnet__Reaction__Element__getName( p_element )
      )
    ); 

  if( p_wd->iZ1 != 0 && p_wd->iA1 != 0 ) {

    p_wd->dScreen *=
      my_pair_screening_function(
        p_wd->dT9,
        p_wd->dRho,
        p_wd->dYe,
        p_wd->iZ1,
        p_wd->iA1,
        i_z2,
        i_a2,
        p_wd->pUserData
      );

  }

  p_wd->iZ1 += i_z2;
  p_wd->iA1 += i_a2;

  return 1;

}

/*##############################################################################
// Base screening function.
//############################################################################*/

double
my_pair_screening_function(
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

  struct user_screening_data user_data = *(struct user_screening_data *) p_data;

  double d_A_bar;

  double d_Gamma_e;
  double d_Gamma_effective;

  double d_H_strong;
  double d_H_weak;
  double d_H_intermediate;

  double d_f = 0.;
  double d_Ye2 = user_data.dYe2;

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
        i_Z1, i_A1, i_Z2, i_A2, d_T9, d_Gamma_e, d_Gamma_effective
      );
    d_H_intermediate = intermediate_screening_factor( d_H_weak, d_H_strong);
    d_f = exp(d_H_intermediate);
  }

  if( (d_Gamma_effective >= 0.8) && (d_Gamma_effective <= 168) )
  {
    d_H_strong =
      strong_screening_factor(
        i_Z1, i_A1, i_Z2, i_A2, d_T9, d_Gamma_e, d_Gamma_effective );
    d_f = exp(d_H_strong);
  }

  if( d_Gamma_effective > 168 ) {
    LIBNUCNET__ERROR(
      "This is beyond the strong screening regime and not covered by this approximation."
    );
  }

  return d_f;

}

/*##############################################################################
// intermediate_screening_factor
//############################################################################*/

double intermediate_screening_factor( double d_H_weak, double d_H_strong ) {

/*=============================================================================
/
/ H_{intermediate} = (H_{weak} H_{strong}) / (H_{weak}^{2} + H_{strong}^2)^{1/2}
/
/=============================================================================*/

  double d_H_intermediate;
  double d_denominator;

  d_denominator = pow(d_H_weak,2.0) + pow(d_H_strong, 2.0);
  d_denominator = pow(d_denominator,1.0/2.0);

  d_H_intermediate =  d_H_weak*d_H_strong/d_denominator;

  return d_H_intermediate;

}


/*##############################################################################
// weak_screening_factor
//############################################################################*/

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

/*=============================================================================
/
/ H_{weak} = Z_{1} Z_{2} \bar{A}^{-1/2} (Y_{e,2} + \theta_{e}Y_{e})^{1/2}
(1.88x10^{8}) (\rho/(T^{3})^{1/2}
/
/ Y_{e,2} = \sum_{i} Z_{i}^{2} Y_{i}
/
/ **WARNINGS**
/
/  1) Equation (A14) in Wallace (1982) contains errors. The above has these errors
corrected.
/  2) The above includes an electron degeneracy factor which would make the
/     Y_{e} term go to zero in degenerate matter. 
/ 
/
/=============================================================================*/

 double d_H_weak;

  d_H_weak = i_Z1*i_Z2*1.88E8;
  d_H_weak = d_H_weak * pow(d_A_bar, -1.0/2.0);
  d_H_weak = d_H_weak * pow(D_THETA_E*d_Ye + d_Ye2, 1.0/2.0);
  d_H_weak = d_H_weak * pow(d_rho, 1.0/2.0) * pow(d_T9*1.0E9, -3.0/2.0);

  return d_H_weak;

}


/*##############################################################################
// strong_screening_factor
//############################################################################*/

double
strong_screening_factor(
  unsigned int i_Z1,
  unsigned int i_A1,
  unsigned int i_Z2,
  unsigned int i_A2,
  double d_T9,
  double d_Gamma_e,
  double d_Gamma_effective
)
{

/*=============================================================================
/
/ \tau_{1,2} = 4.24872((Z_{1}^{2} Z_{2}^{2} \hat{A}_{1,2})/(T_{9}))^{1/3}
/
/ \hat{A}_{1,2} = A_{1}A_{2}\over{A_{1} + A_{2}}
/
/
/ C = 0.896434\Gamma_{effective}\tilde{z} -
/     3.44740(\Gamma_{effective})^{1/4}\tilde{\varsigma} -
/     0.5551(\ln{\Gamma_{e}} + (5/3)\ln{[(Z_{1}Z_{2}/(Z_{1} + Z_{2})]}) - 2.996
/
/ \tilde{z} = (Z_{1} + Z_{2})^{5/3} - Z_{1}^{5/3} - Z_{2}^{5/3}
/
/ \tilde{\varsigma} = (Z_{1} + Z_{2})^{5/12} - Z_{1}^{5/12} - Z_{2}^{5/12}
/
/
/ b = 3.0\Gamma_{effective} / \tau_{1,2}
/
/
/ H_{strong} = C - \tau_{1,2}/3((5/32)b^{3} - 0.014b^{4} -0.128b^{5}) - 
/             \Gamma_{effective}(0.0055b^{4} - 0.0098b^{5} + 0.0048b^{6})
/
/=============================================================================*/

  double d_A_reduced;
  double d_tau12;

  double d_z53;
  double d_z512;

  double d_C;

  double d_b;

  double d_H_strong;

  d_A_reduced = (i_A1*i_A2)/(i_A1 + i_A2);

  d_tau12 = i_Z1*i_Z1*i_Z2*i_Z2*d_A_reduced/d_T9;
  d_tau12 = pow(d_tau12, 1.0/3.0);
  d_tau12 = 4.24872*d_tau12;

  d_z53 =
    pow( (double) (i_Z1 + i_Z2), 5.0/3.0 ) -
    pow( (double) i_Z1, 5.0/3.0) - pow( (double) i_Z2, 5.0/3.0);
  d_z512 = pow( (double) (i_Z1 + i_Z2), 5.0/12.0) -
    pow( (double) i_Z1, 5.0/12.0) - pow( (double) i_Z2, 5.0/12.0);

  d_C = 0.896434 * d_Gamma_e * d_z53;
  d_C -= 3.44740 * pow(d_Gamma_e,1.0/4.0) * d_z512;
  d_C -= 0.5551 * log(d_Gamma_e);
  d_C -= 0.5551 * (5.0/3.0) *
    log( (double) (i_Z1*i_Z2)/ (double) (i_Z1 + i_Z2));
  d_C -= 2.996;

  d_b = 3.0 * d_Gamma_effective / d_tau12;

  d_H_strong =
    d_C - (d_tau12/3.0)*((5.0/32.0)*pow(d_b, 3.0) - 0.014*pow(d_b,4.0) -
    0.128*pow(d_b,5.0));

  d_H_strong += - d_Gamma_effective*(0.0055*pow(d_b,4.0) -
     0.0098*pow(d_b,5.0) + 0.0048*pow(d_b,6.0));

  return d_H_strong; 

}

/*##############################################################################
// calculate_gamma_e
//############################################################################*/

double calculate_gamma_e ( double d_T9, double d_rho, double d_Ye) {

/*=============================================================================
/
/ \Gamma_{e} = (e^2/KT)(4 \pi \rho Y_{e} / 3)^{1/3}
/
/=============================================================================*/

  double d_a_e;
  double d_Gamma_e;

  d_a_e = 4.0*M_PI*d_rho*GSL_CONST_NUM_AVOGADRO*d_Ye/3.0;
  d_a_e = pow(d_a_e,1.0/3.0);

  d_Gamma_e =
    GSL_CONST_CGSM_ELECTRON_CHARGE*GSL_CONST_CGSM_ELECTRON_CHARGE *
    d_a_e / ( GSL_CONST_CGSM_BOLTZMANN * d_T9 * GSL_CONST_NUM_GIGA );

  return d_Gamma_e;

}

/*##############################################################################
// calculate_gamma_effective
//############################################################################*/

double
calculate_gamma_effective(
  unsigned int i_Z1,
  unsigned int i_Z2,
  double d_Gamma_e
)
{

/*=============================================================================
/
/ \Gamma_{effective} = (2 {\over Z_{1} + Z_{2}})^{1/3} Z_{1} Z_{2} \Gamma_{e}
/
/=============================================================================*/

  double d_Gamma_effective;

  d_Gamma_effective = 2.0 / (i_Z1 + i_Z2);
  d_Gamma_effective = pow(d_Gamma_effective,1.0/3.0);
  d_Gamma_effective = d_Gamma_effective * i_Z1 * i_Z2 * d_Gamma_e;

  return d_Gamma_effective;

}
