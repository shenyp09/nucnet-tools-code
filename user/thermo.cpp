////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2015 Clemson University.
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
//! \brief Code for computing relevant thermodynamics quantities.
//!
////////////////////////////////////////////////////////////////////////////////

#include "user/thermo.h"

/**
 * @brief A NucNet Tools namespace for extra (potentially user-supplied)
 *        codes.
 */
namespace user
{

//##############################################################################
// get_thermo_species_list().
//##############################################################################

nnt::species_list_t
get_thermo_species_list( nnt::Zone& zone )
{

  if( zone.hasProperty( nnt::s_THERMO_NUC_VIEW ) )
  {
    return
      nnt::make_species_list(
        Libnucnet__Net__getNuc(
          Libnucnet__NetView__getNet(
            zone.getNetView(
              zone.getProperty<std::string>( nnt::s_THERMO_NUC_VIEW ).c_str()
            )  
          )
        )
      );
  }
  else
  {
    return
      nnt::make_species_list(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      );
  }

}

//############################################################################
// compute_electron_chemical_potential_kT().
//############################################################################

/**
  \brief Compute the electron chemical potential (less the rest mass)
	 divided by kT.

  \param zone A NucNet Tools zone
  \return A double giving the chemical potential / kT.
*/

double
compute_electron_chemical_potential_kT(
  nnt::Zone& zone
)
{

  Libstatmech__Fermion *p_electron;
  double d_result;

  double d_t9 = zone.getProperty<double>( nnt::s_T9 );

  double d_rho = zone.getProperty<double>( nnt::s_RHO );

  double d_ye = Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  d_result =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      d_t9 * GSL_CONST_NUM_GIGA,
      d_rho * d_ye * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_result;

}

//############################################################################
// compute_electron_chemical_potential_kT_temperature_derivative().
//############################################################################

/**
  \brief Compute the temperature derivative of the electron chemical
         potential (less the rest mass) divided by kT.
  \param zone A NucNet Tools zone
  \return A double giving the d(chemical potential / kT)/dT.
*/

double
compute_electron_chemical_potential_kT_temperature_derivative(
  nnt::Zone& zone
)
{


  Libstatmech__Fermion * p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  double d_result =
    Libstatmech__Fermion__computeTemperatureDerivative(
      p_electron,
      S_CHEMICAL_POTENTIAL,
      zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA,
      zone.getProperty<double>( nnt::s_RHO ) *
        GSL_CONST_NUM_AVOGADRO *
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_result;

}

//##############################################################################
// get_electron_neutrino_chemical_potential_kT()
//##############################################################################

/**
 * \brief Returns the neutrino chemical potential for a zone.
 *
 * \return The neutrino chemical potential divided by kT.  If this value is
 *          -infinity, it is returns as GSL_NEGINF.
 */

double
get_electron_neutrino_chemical_potential_kT(
  nnt::Zone& zone
)
{

//  double d_mu_nue_kT;

  try{
    return
      boost::lexical_cast<double>( zone.getProperty<std::string>(
        nnt::s_MU_NUE_KT )
      );
  }
  catch( const boost::bad_lexical_cast& e )
  {
    if( zone.getProperty<std::string>( nnt::s_MU_NUE_KT ) == "-inf" )
    {
      return GSL_NEGINF;
    }
    else
    {
      std::cerr << "Neutrino chemical potential " <<
        zone.getProperty<std::string>( nnt::s_MU_NUE_KT ) <<
        " not valid." << std::endl;
      throw e;
    }
  }

//  return d_mu_nue_kT;

}

//##############################################################################
// compute_baryon_internal_energy_density(). 
//##############################################################################

/**
  \brief Default calculation of the baryon internal energy density.

  \param zone A NucNet Tools zone
  \return A double giving the baryon internal energy density (ergs/cc).
*/

double
compute_baryon_internal_energy_density( nnt::Zone& zone )
{

  double d_eB = 0;
  double d_n, d_T;

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    zone.getProperty<double>( nnt::s_RHO );

  d_T =
    zone.getProperty<double>( nnt::s_T9 ) *
    GSL_CONST_NUM_GIGA;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    d_eB +=
      ( 3. / 2. ) *
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      d_n *
      GSL_CONST_CGSM_BOLTZMANN *
      d_T;
  }

  return d_eB;

}

//##############################################################################
// compute_electron_internal_energy_density().
//##############################################################################

/**
  \brief Default calculation of the electron internal energy density.

  \param zone A NucNet Tools zone
  \return A double giving the electron internal energy density (ergs/cc).
*/

double
compute_electron_internal_energy_density( nnt::Zone& zone )
{

  Libstatmech__Fermion * p_electron;
  double d_T, d_eE = 0;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    nnt::s_INTERNAL_ENERGY_DENSITY,
    0.,
    nnt::d_REL_EPS
  );

  d_T =
    zone.getProperty<double>( nnt::s_T9 ) *
    GSL_CONST_NUM_GIGA;

  d_eE =
    Libstatmech__Fermion__computeQuantity(
      p_electron,
      nnt::s_INTERNAL_ENERGY_DENSITY,
      d_T,
      compute_electron_chemical_potential_kT( zone ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_eE;

}
  
//##############################################################################
// compute_photon_internal_energy_density().
//##############################################################################

/**
  \brief Default calculation of the photon internal energy density.

  \param zone A NucNet Tools zone
  \return A double giving the photon internal energy density (ergs/cc).
*/

double
compute_photon_internal_energy_density( nnt::Zone& zone )
{

  double d_eP = 0, d_T;
  Libstatmech__Boson * p_photon;

  p_photon = Libstatmech__Boson__new( "photon", 0., 2, 0. );

  d_T =
    zone.getProperty<double>( nnt::s_T9 ) *
    GSL_CONST_NUM_GIGA;
  
  d_eP =
    Libstatmech__Boson__computeQuantity(
      p_photon,
      nnt::s_INTERNAL_ENERGY_DENSITY,
      d_T,
      0.,
      NULL,
      NULL
    );

  Libstatmech__Boson__free( p_photon );

  return d_eP;

}

//##############################################################################
// compute_internal_energy_density().
//##############################################################################

/**
  \brief Default calculation of the total internal energy density.

  \param zone A NucNet Tools zone
  \return A double giving the total internal energy density (ergs/cc).
*/

double
compute_internal_energy_density( nnt::Zone& zone )
{

  return
    compute_baryon_internal_energy_density( zone ) +
    compute_electron_internal_energy_density( zone ) +
    compute_photon_internal_energy_density( zone );

}
  
//##############################################################################
// compute_rest_mass_energy_density().
//##############################################################################

/**
 * \brief Compute the total rest mass energy density in a zone.
 * \return The rest mass energy density (cgs).
 */


double
compute_rest_mass_energy_density( nnt::Zone& zone )
{

  double d_n, d_result = 0;

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    zone.getProperty<double>( nnt::s_RHO );

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_result +=
      d_n *
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      (
        Libnucnet__Species__getA( species.getNucnetSpecies() ) *
        GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS *
        gsl_pow_2( GSL_CONST_CGSM_SPEED_OF_LIGHT ) +
        Libnucnet__Species__getMassExcess( species.getNucnetSpecies() ) *
        GSL_CONST_CGSM_ELECTRON_VOLT *
        GSL_CONST_NUM_MEGA
      );
       
  }

  return d_result;

}

//##############################################################################
// compute_baryon_pressure().
//##############################################################################

/**
  \brief Default calculation of the baryon pressure.

  \param zone A NucNet Tools zone
  \return A double giving the baryon pressure (dynes/cm^2).
*/

double
compute_baryon_pressure( nnt::Zone& zone )
{

  double d_pB = 0;
  double d_n, d_T;

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    zone.getProperty<double>( nnt::s_RHO );

  d_T =
    zone.getProperty<double>( nnt::s_T9 ) *
    GSL_CONST_NUM_GIGA;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    d_pB +=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      d_n *
      GSL_CONST_CGSM_BOLTZMANN *
      d_T;
  }

  return d_pB;

}

//##############################################################################
// compute_electron_pressure().
//##############################################################################

/**
  \brief Default calculation of the electron pressure.

  \param zone A NucNet Tools zone
  \return A double giving the electron pressure (dynes/cm^2).
*/

double
compute_electron_pressure( nnt::Zone& zone )
{

  Libstatmech__Fermion * p_electron;
  double d_T, d_pE = 0;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    nnt::s_PRESSURE,
    0.,
    nnt::d_REL_EPS
  );

  d_T = zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA;

  d_pE =
    Libstatmech__Fermion__computeQuantity(
      p_electron,
      nnt::s_PRESSURE,
      d_T,
      compute_electron_chemical_potential_kT( zone ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_pE;

}
  
//##############################################################################
// compute_photon_pressure().
//##############################################################################

/**
  \brief Default calculation of the photon pressure.

  \param zone A NucNet Tools zone
  \return A double giving the photon pressure (dynes/cm^2).
*/

double
compute_photon_pressure( nnt::Zone& zone )
{

  double d_pP = 0, d_T;
  Libstatmech__Boson * p_photon;

  p_photon = Libstatmech__Boson__new( "photon", 0., 2, 0. );

  d_T = zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA;
  
  d_pP =
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_PRESSURE,
      d_T,
      0.,
      NULL,
      NULL
    );

  Libstatmech__Boson__free( p_photon );

  return d_pP;

}

//##############################################################################
// compute_pressure().
//##############################################################################

/**
  \brief Default calculation of the total pressure.

  \param zone A NucNet Tools zone
  \return A double giving the total pressure (dynes/cm^2).
*/

double
compute_pressure( nnt::Zone& zone )
{

  return
    compute_thermo_quantity( zone, nnt::s_PRESSURE, nnt::s_BARYON ) +
    compute_thermo_quantity( zone, nnt::s_PRESSURE, nnt::s_ELECTRON ) +
    compute_thermo_quantity( zone, nnt::s_PRESSURE, nnt::s_PHOTON );

}
  
//##############################################################################
// compute_baryon_dPdT().
//##############################################################################

/**
  \brief Default calculation of the derivative of the baryon pressure with
           temperature.

  \param zone A NucNet Tools zone
  \return A double giving the baryon dPdT (dynes/cm^2/K).
*/

double
compute_baryon_dPdT( nnt::Zone& zone )
{

  double d_B_dPdT = 0;
  double d_n;

  d_n =
    GSL_CONST_NUM_AVOGADRO *
    zone.getProperty<double>( nnt::s_RHO );

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    d_B_dPdT +=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) *
      d_n *
      GSL_CONST_CGSM_BOLTZMANN;
  }

  return d_B_dPdT;

}

//##############################################################################
// compute_electron_dPdT().
//##############################################################################

/**
  \brief Default calculation of the derivative of the electron pressure with
           temperature.

  \param zone A NucNet Tools zone
  \return A double giving the electron dPdT (dynes/cm^2/K).
*/

double
compute_electron_dPdT( nnt::Zone& zone )
{

  Libstatmech__Fermion * p_electron;
  double d_E_dPdT = 0;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    S_PRESSURE,
    0.,
    nnt::d_REL_EPS
  );

  d_E_dPdT =
    Libstatmech__Fermion__computeTemperatureDerivative(
      p_electron,
      S_PRESSURE,
      zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA,
      zone.getProperty<double>( nnt::s_RHO ) *
        GSL_CONST_NUM_AVOGADRO *
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return d_E_dPdT;

}
  
//##############################################################################
// compute_photon_dPdT().
//##############################################################################

/**
  \brief Default calculation of the derivative of the photon pressure with
           temperature.

  \param zone A NucNet Tools zone
  \return A double giving the photon dPdT (dynes/cm^2/K).
*/

double
compute_photon_dPdT( nnt::Zone& zone )
{

  double d_T;

  d_T = zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA;
  
  return
    4. / 3. *
    (
      4. *
      GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT /
      GSL_CONST_CGSM_SPEED_OF_LIGHT
    ) *
    gsl_pow_3( d_T );

}

//##############################################################################
// compute_dPdT().
//##############################################################################

/**
  \brief Default calculation of the derivative of the total pressure with
           temperature.

  \param zone A NucNet Tools zone
  \return A double giving the total dPdT (dynes/cm^2/K).
*/

double
compute_dPdT( nnt::Zone& zone )
{

  return
    compute_baryon_dPdT( zone ) +
    compute_electron_dPdT( zone ) +
    compute_photon_dPdT( zone );

}
  
//##############################################################################
// compute_baryon_entropy_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the baryon entropy per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the baryon entropy (per k_B per nucleon).
*/

double
compute_baryon_entropy_per_nucleon( nnt::Zone& zone )
{

  double d_result = 0;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    if( 
	Libnucnet__Zone__getSpeciesAbundance(
	  zone.getNucnetZone(),
	  species.getNucnetSpecies()
	) > 1.e-30
    )
    {

      d_result +=
	(
	  2.5 *
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    species.getNucnetSpecies()
	  )
	)
	-
	(
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    species.getNucnetSpecies()
	  )
	  *
	  log(
	    Libnucnet__Zone__getSpeciesAbundance(
	      zone.getNucnetZone(),
	      species.getNucnetSpecies()
	    ) /
	    Libnucnet__Species__computeQuantumAbundance(
	      species.getNucnetSpecies(),
	      zone.getProperty<double>( nnt::s_T9 ),
	      zone.getProperty<double>( nnt::s_RHO )
	    )
	  )
	);

    }

  }

  return d_result;

}

//##############################################################################
// compute_electron_entropy_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the electron entropy per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the electron entropy (per k_B per nucleon).
*/

double
compute_electron_entropy_per_nucleon( nnt::Zone& zone )
{

  double d_result = 0;
  Libstatmech__Fermion * p_electron;
  double d_T;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    S_ENTROPY_DENSITY,
    0.,
    nnt::d_REL_EPS
  );

  d_T = zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA;

  d_result =
    Libstatmech__Fermion__computeQuantity(
      p_electron,
      S_ENTROPY_DENSITY,
      d_T,
      compute_electron_chemical_potential_kT( zone ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return
    d_result /
    (
      zone.getProperty<double>( nnt::s_RHO ) *
      GSL_CONST_NUM_AVOGADRO *
      GSL_CONST_CGSM_BOLTZMANN
    );

}
  
//##############################################################################
// compute_photon_entropy_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the photon entropy per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the photon entropy (per k_B per nucleon).
*/

double
compute_photon_entropy_per_nucleon( nnt::Zone& zone )
{

  double d_T, d_result = 0;
  Libstatmech__Boson * p_photon;

  p_photon = Libstatmech__Boson__new( "photon", 0., 2, 0. );

  d_T = zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA;
  
  d_result =
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_ENTROPY_DENSITY,
      d_T,
      0.,
      NULL,
      NULL
    );

  Libstatmech__Boson__free( p_photon );

  return
    d_result /
    (
      zone.getProperty<double>( nnt::s_RHO ) *
      GSL_CONST_NUM_AVOGADRO *
      GSL_CONST_CGSM_BOLTZMANN
    );


}

//##############################################################################
// compute_entropy_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the total entropy per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the total entropy (per k_B per nucleon).
*/

double
compute_entropy_per_nucleon( nnt::Zone& zone )
{

  return
    compute_baryon_entropy_per_nucleon( zone ) +
    compute_electron_entropy_per_nucleon( zone ) +
    compute_photon_entropy_per_nucleon( zone );

}
  
//##############################################################################
// compute_baryon_specific_heat_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the baryon specific heat per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the baryon specific heat (per k_B per K per nucleon).
*/

double
compute_baryon_specific_heat_per_nucleon( nnt::Zone& zone )
{

  double d_cvb = 0;

  nnt::species_list_t species_list = get_thermo_species_list( zone );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    d_cvb +=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      )
      *
      (
        1.5
        +
        (
          zone.getProperty<double>( nnt::s_T9 ) *
          GSL_CONST_NUM_GIGA *
          compute_dlnG_dT(
            species.getNucnetSpecies(), 
            zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA
          )
        )
      );

  }

  return d_cvb;

}

//##############################################################################
// compute_photon_specific_heat_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the photon specific heat per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the photon specific heat (per k_B per K per nucleon).
*/

double
compute_photon_specific_heat_per_nucleon( nnt::Zone& zone )
{

  return
    4.
    *
    (
      4. * GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT /
      GSL_CONST_CGSM_SPEED_OF_LIGHT
    )
    *
    gsl_pow_3( zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA )
    /
    (
      zone.getProperty<double>( nnt::s_RHO ) *
        GSL_CONST_NUM_AVOGADRO *
        GSL_CONST_CGSM_BOLTZMANN
    );

}

//##############################################################################
// compute_electron_specific_heat_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the electron specific heat per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the electron specific heat
            (per k_B per K per nucleon).
*/

double
compute_electron_specific_heat_per_nucleon( nnt::Zone& zone )
{
  
  Libstatmech__Fermion * p_electron;
  double d_T, d_cve = 0, d_ne;

  p_electron = 
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON, nnt::d_ELECTRON_MASS_IN_MEV, 2, -1.
    );

  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_electron,
    S_ENTROPY_DENSITY,
    0.,
    nnt::d_REL_EPS
  );

  d_T = zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA;

  d_ne =
    zone.getProperty<double>( nnt::s_RHO ) *
    GSL_CONST_NUM_AVOGADRO *
    Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

  d_cve =
    d_T *
    Libstatmech__Fermion__computeTemperatureDerivative(
      p_electron,
      S_ENTROPY_DENSITY,
      d_T,
      d_ne,
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  return
    d_cve /
    (
      zone.getProperty<double>( nnt::s_RHO ) *
      GSL_CONST_NUM_AVOGADRO *
      GSL_CONST_CGSM_BOLTZMANN
    );

}
  
//##############################################################################
// compute_specific_heat_per_nucleon().
//##############################################################################

/**
  \brief Default calculation of the total specific heat per nucleon.

  \param zone A NucNet Tools zone
  \return A double giving the total specific heat (per k_B per K per nucleon).
*/

double
compute_specific_heat_per_nucleon( nnt::Zone& zone )
{

  return
    compute_baryon_specific_heat_per_nucleon( zone ) +
    compute_electron_specific_heat_per_nucleon( zone ) +
    compute_photon_specific_heat_per_nucleon( zone );

}
  
//##############################################################################
// assign_default_thermo_functions().
//##############################################################################

void
assign_default_thermo_functions( nnt::Zone& zone )
{

  //============================================================================
  // Pressure.
  //============================================================================
    
  if( !zone.hasFunction( nnt::char_cat( nnt::s_BARYON, nnt::s_PRESSURE ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_PRESSURE ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_baryon_pressure
      ),
      "The baryon pressure in dynes / cm^2.",
      nnt::s_THERMO
    );
  }

  if( !zone.hasFunction( nnt::char_cat( nnt::s_ELECTRON, nnt::s_PRESSURE ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_PRESSURE ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_electron_pressure
      ),
      "The electron pressure in dynes / cm^2.",
      nnt::s_THERMO
    );
  }

  if( !zone.hasFunction( nnt::char_cat( nnt::s_PHOTON, nnt::s_PRESSURE ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_PRESSURE ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_photon_pressure
      ),
      "The photon pressure in dynes / cm^2.",
      nnt::s_THERMO
    );
  }

  if( !zone.hasFunction( nnt::char_cat( nnt::s_TOTAL, nnt::s_PRESSURE ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_PRESSURE ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_pressure
      ),
      "The total pressure in dynes / cm^2.",
      nnt::s_THERMO
    );
  }

  //============================================================================
  // Entropy per nucleon.
  //===========================================================================
  
  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_ENTROPY_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_ENTROPY_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_baryon_entropy_per_nucleon
      ),
      "The baryon entropy per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_ENTROPY_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_ENTROPY_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_electron_entropy_per_nucleon
      ),
      "The electron entropy per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_ENTROPY_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_ENTROPY_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_photon_entropy_per_nucleon
      ),
      "The photon entropy per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_ENTROPY_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_ENTROPY_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_entropy_per_nucleon
      ),
      "The total entropy per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  //============================================================================
  // dPdT.
  //============================================================================
  
  if( !zone.hasFunction( nnt::char_cat( nnt::s_BARYON, nnt::s_DPDT ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_DPDT ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_baryon_dPdT
      ),
      "The derivative of baryon pressure with respect to temperature (in dynes per cm^2 per Kelvin",
      nnt::s_THERMO
    );
  }

  if( !zone.hasFunction( nnt::char_cat( nnt::s_ELECTRON, nnt::s_DPDT ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_DPDT ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_electron_dPdT
      ),
      "The derivative of electron pressure with respect to temperature (in dynes per cm^2 per Kelvin",
      nnt::s_THERMO
    );
  }

  if( !zone.hasFunction( nnt::char_cat( nnt::s_PHOTON, nnt::s_DPDT ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_DPDT ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_photon_dPdT
      ),
      "The derivative of photon pressure with respect to temperature (in dynes per cm^2 per Kelvin",
      nnt::s_THERMO
    );
  }

  if( !zone.hasFunction( nnt::char_cat( nnt::s_TOTAL, nnt::s_DPDT ) ) )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_DPDT ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_dPdT
      ),
      "The derivative of pressure with respect to temperature (in dynes per cm^2 per Kelvin",
      nnt::s_THERMO
    );
  }

  //============================================================================
  // Specific heat.
  //============================================================================
  
  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_SPECIFIC_HEAT_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_SPECIFIC_HEAT_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_baryon_specific_heat_per_nucleon
      ),
      "The baryon specific heat per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_SPECIFIC_HEAT_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_SPECIFIC_HEAT_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_electron_specific_heat_per_nucleon
      ),
      "The electron specific heat per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_SPECIFIC_HEAT_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_SPECIFIC_HEAT_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_photon_specific_heat_per_nucleon
      ),
        "The photon specific heat per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_SPECIFIC_HEAT_PER_NUCLEON )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_SPECIFIC_HEAT_PER_NUCLEON ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_specific_heat_per_nucleon
      ),
      "The total specific heat per nucleon in units of Boltzmann's constant.",
      nnt::s_THERMO
    );
  }

  //============================================================================
  // Internal energy density.
  //============================================================================
  
  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_INTERNAL_ENERGY_DENSITY )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_BARYON, nnt::s_INTERNAL_ENERGY_DENSITY ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_baryon_internal_energy_density
      ),
      "The baryon internal energy density in ergs per cm^3.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_INTERNAL_ENERGY_DENSITY )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_INTERNAL_ENERGY_DENSITY ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_electron_internal_energy_density
      ),
      "The electron internal energy density in ergs per cm^3.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_INTERNAL_ENERGY_DENSITY )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_PHOTON, nnt::s_INTERNAL_ENERGY_DENSITY ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_photon_internal_energy_density
      ),
      "The photon internal energy density in ergs per cm^3.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_INTERNAL_ENERGY_DENSITY )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_TOTAL, nnt::s_INTERNAL_ENERGY_DENSITY ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_internal_energy_density
      ),
      "The total internal energy density in ergs per cm^3.",
      nnt::s_THERMO
    );
  }

  //============================================================================
  // Chemical potential / kT.
  //============================================================================
  
  if(
    !zone.hasFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_CHEMICAL_POTENTIAL_KT )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_ELECTRON, nnt::s_CHEMICAL_POTENTIAL_KT ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_electron_chemical_potential_kT
      ),
      "The electron chemical potential (less rest mass) divided by kT.",
      nnt::s_THERMO
    );
  }

  if(
    !zone.hasFunction(
      nnt::char_cat(
        nnt::s_ELECTRON,
        nnt::s_T_DERIVATIVE_CHEMICAL_POTENTIAL_KT
      )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat(
        nnt::s_ELECTRON,
        nnt::s_T_DERIVATIVE_CHEMICAL_POTENTIAL_KT
      ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        compute_electron_chemical_potential_kT_temperature_derivative
      ),
      "The temperature derivative of the electron chemical potential (less rest mass) divided by kT (in units of K^-1).",
      nnt::s_THERMO
    );
  }

  if(
     !zone.hasFunction(
       nnt::char_cat( nnt::s_NEUTRINO_E, nnt::s_CHEMICAL_POTENTIAL_KT )
    )
  )
  {
    zone.updateFunction(
      nnt::char_cat( nnt::s_NEUTRINO_E, nnt::s_CHEMICAL_POTENTIAL_KT ),
      static_cast<boost::function<double( nnt::Zone& )> >(
        get_electron_neutrino_chemical_potential_kT
      ),
      "The electron neutrino chemical potential divided by kT.",
      nnt::s_THERMO
    );
  }

}

//##############################################################################
// list_zone_thermo_quantities().
//##############################################################################

/**
 * \brief List the thermo quantities available to compute for a zone.
 * \param zone A NucNet zone for which the quantities are defined.
 * \return A vector of strings the currently available quantities.
 */

std::vector<std::string>
list_zone_thermo_quantities( nnt::Zone& zone )
{

  assign_default_thermo_functions( zone );

  std::vector<std::string> quantities =
    zone.getListOfFunctions( nnt::s_THERMO );

  std::sort( quantities.begin(), quantities.end() );

  return quantities;

}
  
//##############################################################################
// get_zone_thermo_quantity_doc().
//##############################################################################

/**
 * \brief Get the documentation of a thermo quantity for a zone.
 * \param zone The zone for which the quantity is defined.
 * \param s_quantity A string giving the quantity.
 * \return A string of the documentation.
 */

std::string
get_zone_thermo_quantity_doc(
  nnt::Zone& zone,
  const std::string& s_quantity
)
{

  if( !zone.hasFunction( s_quantity ) )
  {
    assign_default_thermo_functions( zone );
  }

  return zone.getFunctionDoc( s_quantity );

}
  
//##############################################################################
// compute_thermo_quantity().
//##############################################################################

/**
 * \brief Compute a thermodynamic quantity from a zone.
 * \param zone A Nucnet-Tools zone.
 * \param s_quantity A string giving the quantity to compute.
 * \return The computed quantity.
 */

double
compute_thermo_quantity(
  nnt::Zone& zone,
  const std::string& s_quantity
)
{

  if( !zone.hasFunction( s_quantity ) )
  {
    assign_default_thermo_functions( zone );
  }

  if( zone.hasFunction( s_quantity ) )
  {
    return
      boost::any_cast<boost::function<double( nnt::Zone& )> >(
        zone.getFunction( s_quantity )
      )( zone );
  }
  else
  {
    std::cerr <<
      std::endl << s_quantity << " not valid" << std::endl;
    exit( EXIT_FAILURE );
  }

}
  
//##############################################################################
// compute_thermo_quantity().
//##############################################################################

/**
 * \brief Compute a thermodynamic quantity from a zone.
 * \param zone A Nucnet-Tools zone.
 * \param s_quantity A string giving the quantity to compute.
 * \param s_particle The particle ("baryon", "electron", or "photon")
 *     to compute, or the total ("total").
 * \return The computed quantity.
 */

double
compute_thermo_quantity(
  nnt::Zone& zone,
  const std::string& s_quantity,
  const std::string& s_particle
)
{

  return compute_thermo_quantity( zone, s_particle + " " + s_quantity );

}
  
//##############################################################################
// compute_log10_t9_entropy_root_with_equilibrium(). 
//##############################################################################

double
compute_log10_t9_entropy_root_with_equilibrium(
  double d_log10_t9,
  nnt::Zone& zone
)
{

  zone.updateProperty(
    nnt::s_T9,
    pow( 10., d_log10_t9 )
  );

  nnt::set_zone_abundances_to_equilibrium( zone );

  return
    compute_thermo_quantity(
      zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    ) -
    zone.getProperty<double>( nnt::s_ENTROPY_PER_NUCLEON );

}
    
//##############################################################################
// compute_log10_t9_entropy_root().
//##############################################################################

double
compute_log10_t9_entropy_root(
  double d_log10_t9,
  nnt::Zone& zone
)
{

  zone.updateProperty(
    nnt::s_T9,
    pow( 10., d_log10_t9 )
  );

  return
    compute_thermo_quantity(
      zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    ) -
    zone.getProperty<double>( nnt::s_ENTROPY_PER_NUCLEON );

}

//##############################################################################
// compute_log10_density_entropy_root_with_equilibrium(). 
//##############################################################################

double
compute_log10_density_entropy_root_with_equilibrium(
  double d_log10_rho,
  nnt::Zone& zone 
)
{

  zone.updateProperty(
    nnt::s_RHO,
    pow( 10., d_log10_rho )
  );

  nnt::set_zone_abundances_to_equilibrium( zone );

  return
    compute_thermo_quantity(
      zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    )
    -
    zone.getProperty<double>( nnt::s_ENTROPY_PER_NUCLEON );

}

//##############################################################################
// compute_specific_heat_per_nucleon_at_constant_pressure().
//##############################################################################

double
compute_specific_heat_per_nucleon_at_constant_pressure( nnt::Zone& zone )
{

  zone.updateProperty(
    nnt::s_PRESSURE,
    compute_thermo_quantity(
      zone,
      nnt::s_PRESSURE,
      nnt::s_TOTAL
    )
  );

  double d_t9 = zone.getProperty<double>( nnt::s_T9 );
  
  double d_cp =
    nnt::compute_derivative(
      boost::bind(
        dS_dT_at_constant_pressure,
        _1,
        boost::ref( zone )
      ),
      zone.getProperty<double>( nnt::s_T9 )
    );

  zone.updateProperty( nnt::s_T9, d_t9 );

  return d_t9 * d_cp;

}
        
//##############################################################################
// dS_dT_at_constant_pressure().
//##############################################################################

double dS_dT_at_constant_pressure(
  double d_T9,
  nnt::Zone& zone
)
{

  zone.updateProperty(
    nnt::s_T9, 
    d_T9
  );

  zone.updateProperty(
    nnt::s_RHO, 
    pow(
      10.,
      nnt::compute_1d_root(
        boost::bind( compute_log10_rho_pressure_root, _1, boost::ref( zone ) ),
        log10( zone.getProperty<double>( nnt::s_RHO ) ),
        2
      )
    )
  );

  return
    compute_thermo_quantity(
      zone,
      nnt::s_TOTAL,
      nnt::s_ENTROPY_PER_NUCLEON
    );

}

//##############################################################################
// compute_log10_rho_pressure_root().
//##############################################################################

double
compute_log10_rho_pressure_root(
  double d_log10_rho,
  nnt::Zone& zone
)
{

  zone.updateProperty(
    nnt::s_RHO,
    pow( 10., d_log10_rho )
  );

  double d_result =
    compute_pressure( zone )
    -
    zone.getProperty<double>( nnt::s_PRESSURE );

  return d_result;

}

//##############################################################################
// compute_sound_speed().
//##############################################################################

double
compute_sound_speed( nnt::Zone& zone )
{

  double d_cs;

  zone.updateProperty(
    nnt::s_ENTROPY_PER_NUCLEON,
    compute_thermo_quantity(
    zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    )
  );

  std::string s_rho = zone.getProperty<std::string>( nnt::s_RHO );
  
  d_cs =
    nnt::compute_derivative(
      boost::bind(
        dP_drho,
        _1,
        boost::ref( zone )
      ),
      zone.getProperty<double>( nnt::s_RHO )
    );

  zone.updateProperty( nnt::s_RHO, s_rho );

  return sqrt( d_cs );

}
        
//##############################################################################
// dP_drho().
//##############################################################################

double dP_drho(
  double d_rho,
  nnt::Zone& zone
)
{

  double d_t9;

  std::string s_t9 = zone.getProperty<std::string>( nnt::s_T9 );

  zone.updateProperty(
    nnt::s_RHO, 
    d_rho
  );

  d_t9 =
    pow(
      10.,
      nnt::compute_1d_root(
        boost::bind( compute_log10_t9_entropy_root, _1, boost::ref( zone ) ),
        log10( zone.getProperty<double>( nnt::s_T9 ) ),
        1.1
      )
    );

  zone.updateProperty(
    nnt::s_T9,
    d_t9 
  );

  double d_result = compute_pressure( zone );

  zone.updateProperty( nnt::s_T9, s_t9 );

  return d_result;

}

//##############################################################################
// compute_dlnG_dT()
//##############################################################################

double
compute_dlnG_dT(
  Libnucnet__Species * p_species,
  double d_T
)
{

  return
    nnt::compute_derivative(
      boost::bind(
        compute_lnG,
        _1,
        p_species
      ),
      d_T
    );

}

//##############################################################################
// compute_lnG().
//##############################################################################

double
compute_lnG( double d_T, Libnucnet__Species * p_species )
{

  return
    log(
      Libnucnet__Species__computePartitionFunction(
        p_species,
        d_T / GSL_CONST_NUM_GIGA
      )
    );

}

} // namespace user
