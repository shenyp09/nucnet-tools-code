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
//! \brief Code for user-defined NSE correction functions.
//!
////////////////////////////////////////////////////////////////////////////////

#include "user/nse_corr.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */
namespace user
{

//##############################################################################
//   In NucNet Tools, the NSE correction factor corrects for deviations
// away from the ideal gas expression for the chemical potential.  The
// ideal, non-relativistic expression for the chemical potential mu_i of a
// gas of species i is
//
//   \mu_i = m_ic^2 + kT \ln\left( \frac{Y_i}{Y_{Qi}} \right)                (1)
//
// where m_i is the rest mass of species i, k is Boltzmann's constant,
// T is the temperature, Y_i is the abundance per nucleon of species i,
// and Y_{Qi} is the quantum abundance of species i.  The quantum abundance
// is given by
//
//   Y_{Qi} =
//     \frac{G_i}{\rho N_A} \left( \frac{m_i kT}{2\pi\hbar^2} \right)^{3/2}  (2)
//
// where G_i is the partition function for species i, \rho is the mass density,
// and N_A is Avogadro's number.  The chemical potential less the rest mass
// mu_i' for the ideal, non-relativistic gas case is then
//
//   \mu_i' = kT \ln\left( \frac{Y_i}{Y_{Qi}} \right)                        (3)
//
// The NSE correction factor f_{corr,i} for species i is defined as
//
//   f_{corr,i} =
//     \ln\left( \frac{Y_i}{Y_{Qi}} \right) - \frac{\mu_i'}{kT}              (4)
//
// If the chemical potential for the system at hand is
//
//   \mu_i' = kT \ln\left( \frac{Y_i}{Y_{Qi}} \right) + \Delta \mu_i         (5)
//
// where \Delta mu_i is the deviation from the ideal gas expression, then from
// (5), it is then evident that
//
//   f_{corr,i} = -\frac{\Delta \mu_i}{kT}                                   (6)
//
// For more details, the user should consult the libnucnet technical report
// "Screening and Reverse Rate Correction Factors in libnucnet" available
// at http://libnucnet.sourceforge.net.
//
// In this example file, we consider Coulomb corrections to the NSE.  The
// correction to the chemical potential is \mu_{i,C}, which accounts for
// the interaction of the species with the surrounding plasma.  The NSE
// corection factor is thus
// 
//   f_{corr,i} = -\frac{\mu_{i,C}}{kT}                                      (7)
//
// The Coulomb correction used is that of Bravo and Garcia-Senz (MNRAS,
// 307, 984-992, (1999)).  The Helmholtz free energy / kT per particle
// is their Eq. (1)
//
//   \frac{f_C}(Gamma_i)}{kT} =
//    a Gamma_i + 4b Gamma_i^{1/4} - 4c Gamma_i^{-1/4} + d ln( Gamma_i ) - o (8)
//
// with
//
//   o = a + 4(b - c) - \frac{f_C(Gamma_i = 1)}{kT}                          (9)
//
// and where
//
//   Gamma_i = Z^{5/3} Gamma_e                                              (10)
//
// with Z the species i charge and Gamma_e defined below.  The constants
// a, b, c, and d are from Ogata and Ichimaru (Phys. Rev. A 36, 5451, (1987)).
// We use the Gamma_i < 1 expression suggested by Bravo and Garcia-Senz
// and match the value and the first derivative of the Helmholtz free
// energy per ion for Gamma_i < 1 and Gamma_i > 1 expressions at Gamma_i = 1.
//##############################################################################

//##############################################################################
// set_nse_correction_function().
//##############################################################################

void
set_nse_correction_function( nnt::Zone& zone )
{

  Libnucnet__Zone__setNseCorrectionFactorFunction(
    zone.getNucnetZone(),
    (Libnucnet__Species__nseCorrectionFactorFunction) nse_correction,
    NULL
  );

  zone.updateFunction(
    nnt::s_NSE_CORRECTION_FACTOR_DATA_FUNCTION,
    static_cast<boost::function<boost::any()> >(
      boost::bind( get_nse_correction_data )
    )
  );
     
}

//##############################################################################
// Routine to get the user-supplied Coulomb correction factor data.
//##############################################################################

nse_corr_data_t
get_nse_correction_data( )
{

  nse_corr_data_t corr_data;

  return corr_data;

}

//##############################################################################
// Gamma_e
//##############################################################################

double Gamma_e ( double d_T9, double d_rho, double d_Ye)
{

//=============================================================================
//
// \Gamma_{e} = (e^2/a_e kT)
//
// a_e = ( 4 * \pi * \rho * N_A * Y_e / 3 )^{1/3}
//
// a_e is the electron cloud radius.  Notice that gsl gives the electron
// charge e in abampere-s.  To convert to esu, multiply by the speed of light.
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
// Species Coulomb chemical potential.
//##############################################################################

double
species_coulomb_chemical_potential(
  Libnucnet__Species *p_species,
  double d_t9,
  double d_rho,
  double d_ye,
  void * p_data
)
{

  unsigned int i_z;
  double Gamma_i;

  if( !p_species )
  {
    std::cerr << "Invalid species!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !p_data )
  {
    std::cerr << "Missing extra data for nse_correction()." << std::endl;
    exit( EXIT_FAILURE );
  }

  nse_corr_data_t nse_corr_data =
    boost::any_cast<nse_corr_data_t>( *(boost::any *) p_data );

  i_z = Libnucnet__Species__getZ( p_species );

  if( i_z == 0 ) { return 0; }   // No correction for the neutron.

  Gamma_i = Gamma_e( d_t9, d_rho, d_ye ) * pow( (double) i_z, 5. / 3. );

  if( Gamma_i > 1 )
  {
    return
      nse_corr_data.a * Gamma_i
      +
      4. * nse_corr_data.b * pow( Gamma_i, 1. / 4. )
      -
      4. * nse_corr_data.c * pow( Gamma_i, -1. / 4. )
      +
      nse_corr_data.d * log( Gamma_i )
      -
      nse_corr_data.o;
  }
  else
  {
    return
      nse_corr_data.beta * pow( Gamma_i, nse_corr_data.gamma ) /
        nse_corr_data.gamma
      -
      pow( Gamma_i, 3. / 2. ) / sqrt(3.);
  }

}

//##############################################################################
// User-supplied Coulomb correction factor function based on Eq. (12) of
// Bravo and Garcia-Senz (1999).  It is the negative of muiC/kT.
//##############################################################################

double
nse_correction(
  Libnucnet__Species *p_species,
  double d_t9,
  double d_rho,
  double d_ye,
  void * p_data
)
{

  return
    -species_coulomb_chemical_potential( p_species, d_t9, d_rho, d_ye, p_data );

}

//##############################################################################
// Species Coulomb energy per particle (in units of kT).
//##############################################################################

double
species_coulomb_energy(
  Libnucnet__Species *p_species,
  double d_t9,
  double d_rho,
  double d_ye,
  void * p_data
)
{

  unsigned int i_z;
  double Gamma_i;

  if( !p_species )
  {
    std::cerr << "Invalid species!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !p_data )
  {
    std::cerr << "Missing extra data for nse_correction()." << std::endl;
    exit( EXIT_FAILURE );
  }

  nse_corr_data_t nse_corr_data =
    boost::any_cast<nse_corr_data_t>( *(boost::any *) p_data );

  i_z = Libnucnet__Species__getZ( p_species );

  if( i_z == 0 ) { return 0; }   // No correction for the neutron.

  Gamma_i = Gamma_e( d_t9, d_rho, d_ye ) * pow( (double) i_z, 5. / 3. );

  if( Gamma_i > 1 )
  {
    return
      nse_corr_data.a * Gamma_i
      +
      nse_corr_data.b * pow( Gamma_i, 1. / 4. )
      +
      nse_corr_data.c * pow( Gamma_i, -1. / 4. )
      +
      nse_corr_data.d;
  }
  else
  {
    return
      nse_corr_data.beta * pow( Gamma_i, nse_corr_data.gamma )
      -
      sqrt( 3. ) * pow( Gamma_i, 3. / 2. ) / 2.;
  }

}

//##############################################################################
// Species Coulomb entropy per particle (in units of Boltzmann's constant).
//##############################################################################

double
species_coulomb_entropy(
  Libnucnet__Species *p_species,
  double d_t9,
  double d_rho,
  double d_ye,
  void * p_data
)
{

  unsigned int i_z;
  double Gamma_i;

  if( !p_species )
  {
    std::cerr << "Invalid species!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( !p_data )
  {
    std::cerr << "Missing extra data for nse_correction()." << std::endl;
    exit( EXIT_FAILURE );
  }

  nse_corr_data_t nse_corr_data =
    boost::any_cast<nse_corr_data_t>( *(boost::any *) p_data );

  i_z = Libnucnet__Species__getZ( p_species );

  if( i_z == 0 ) { return 0; }   // No correction for the neutron.

  Gamma_i = Gamma_e( d_t9, d_rho, d_ye ) * pow( (double) i_z, 5. / 3. );

  if( Gamma_i > 1 )
  {
    return
     -3. * nse_corr_data.b * pow( Gamma_i, 1. / 4. )
      +
      5. * nse_corr_data.c * pow( Gamma_i, -1. / 4. )
      +
      nse_corr_data.d * ( 1. - log( Gamma_i ) )
      +
      nse_corr_data.o;
  }
  else
  {
    return
      nse_corr_data.beta * pow( Gamma_i, nse_corr_data.gamma ) *
        ( 1. - 1. / nse_corr_data.gamma )
      -
      pow( Gamma_i, 3. / 2. ) / ( 2. * sqrt(3.) );
  }

}

} // namespace user
