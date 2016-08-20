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
//! \brief Example code to compute the temperature from a given density
//!    and entropy per nucleon in an equilibrium.
//!
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Include's.
//##############################################################################

#include <Libstatmech.h>
#include <Libnucnet.h>
#include <Libnuceq.h>

#include "nnt/string_defs.h"
#include "nnt/auxiliary.h"

#include "user/thermo.h"

//##############################################################################
// Defines.
//##############################################################################

#define MY_BUF_SIZE  32

//##############################################################################
// count_clusters().
//##############################################################################

void
count_clusters(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value,
  size_t * p_clusters
)
{

  if( s_name && s_tag1 && !s_tag2 && s_value ) { (*p_clusters)++; }

}

//##############################################################################
// insert_data().
//##############################################################################

void
insert_data(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value,
  std::vector<std::string> * p_v
)
{

  if( s_name && s_tag1 && !s_tag2 && s_value )
  { 
    (*p_v)[boost::lexical_cast<size_t>( s_tag1 )] = s_value;
  }

}

//##############################################################################
// print_equilibrium().
//##############################################################################

void
print_equilibrium( nnt::Zone& wrapped_zone )
{

  size_t i_clusters = 0;

  Libnucnet__Zone__iterateOptionalProperties(
    wrapped_zone.getNucnetZone(),
    "cluster xpath",
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) count_clusters,
    &i_clusters
  );

  std::vector<std::string> v_xpath( i_clusters ), v_mukT( i_clusters );

  Libnucnet__Zone__iterateOptionalProperties(
    wrapped_zone.getNucnetZone(),
    nnt::s_CLUSTER_XPATH,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) insert_data,
    &v_xpath
  );
    
  Libnucnet__Zone__iterateOptionalProperties(
    wrapped_zone.getNucnetZone(),
    nnt::s_CLUSTER_CHEMICAL_POTENTIAL,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) insert_data,
    &v_mukT
  );

  for( size_t i = 0; i < i_clusters; i++ )
  {
    std::cout << "Cluster: " << v_xpath[i] << "   " <<
      "mu_cluster / kT: " << v_mukT[i] << std::endl;
  }
    
  std::cout << std::endl;

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;
  nnt::Zone wrapped_zone;
  int i_count;
  char s_property[MY_BUF_SIZE];

  //============================================================================
  // Check input.
  //============================================================================

  if( argc < 5 )
  {
    fprintf(
      stderr,
      "Usage: argv[0] input_net_xml rho s/b Ye [xpath constraint]\n\n"
    );
    fprintf(
      stderr,
      "  input_net_xml = input net xml file\n\n"
    );
    fprintf(
      stderr,
      "  rho = input mass density (g/cc)\n\n"
    );
    fprintf(
      stderr,
      "  s/b = input entropy per nucleon\n\n"
    );
    fprintf(
      stderr,
      "  out_xml = output xml file\n\n"
    );
    fprintf(
      stderr,
      "  Ye = input Ye (optional)\n\n"
    );
    fprintf(
      stderr,
      "  [xpath constraint] = xpath for equilibrium cluster and abundance constraint (optional)\n\n"
    );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Get input.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  p_zone =
    Libnucnet__Zone__new( Libnucnet__getNet( p_my_nucnet ), "0", "0", "0" );

  Libnucnet__addZone( p_my_nucnet, p_zone );

  //============================================================================
  // Input.
  //============================================================================

  Libnucnet__Zone__updateProperty( p_zone, nnt::s_RHO, NULL, NULL, argv[2] );

  Libnucnet__Zone__updateProperty(
    p_zone, 
    nnt::s_ENTROPY_PER_NUCLEON,
    NULL,
    NULL,
    argv[3]
  );

  //============================================================================
  // Optional input.
  //============================================================================

  if( argc > 5 )
    Libnucnet__Zone__updateProperty( p_zone, nnt::s_YE, NULL, NULL, argv[5] );

  if( argc > 6 )
  {
    i_count = 6;
    while( i_count < argc )
    {

      sprintf( s_property, "%d", ( i_count - 6 ) / 2 );
      Libnucnet__Zone__updateProperty(
        p_zone,
        nnt::s_CLUSTER_XPATH,
        s_property,
        NULL,
        argv[i_count++]
      );
      Libnucnet__Zone__updateProperty(
        p_zone,
        nnt::s_CLUSTER_CONSTRAINT,
        s_property,
        NULL,
        argv[i_count++]
      );
    }

  }

  //============================================================================
  // Set wrapped zone.
  //============================================================================

  wrapped_zone.setNucnetZone( p_zone );

  //============================================================================
  // Compute t9 from entropy.
  //============================================================================

  wrapped_zone.updateProperty(
    nnt::s_T9,
    pow(
      10.,
      nnt::compute_1d_root(
        boost::bind(
          user::compute_log10_t9_entropy_root_with_equilibrium,
          _1,
          boost::ref( wrapped_zone )
        ),
        1.,
        1.1
      )
    )
  );
      
  //============================================================================
  // Print out basic data.
  //============================================================================

  fprintf(
    stdout,
    "\nTemperature (K) = %g\n",
    wrapped_zone.getProperty<double>( nnt::s_T9 ) * GSL_CONST_NUM_GIGA
  );

  fprintf(
    stdout,
    "Density = %s g/cc, Ye = %g\n\n",
    argv[2],
    Libnucnet__Zone__computeZMoment( p_zone, 1 )
  );

  //============================================================================
  // Print out equilibrium cluster data.
  //============================================================================

  print_equilibrium( wrapped_zone );
     
  //============================================================================
  // Print out entropy components.
  //============================================================================

  std::cout <<
    "Baryon entropy per nucleon = " <<
    user::compute_thermo_quantity(
      wrapped_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_BARYON
    ) <<
    std::endl;
     
  std::cout << std::endl <<
    "Electron entropy per nucleon = " <<
    user::compute_thermo_quantity(
      wrapped_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_ELECTRON
    ) <<
    std::endl;
     
  std::cout << std::endl <<
    "Photon entropy per nucleon = " <<
    user::compute_thermo_quantity(
      wrapped_zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_PHOTON
    ) <<
    std::endl;

  std::cout << std::endl;
     
  //============================================================================
  // Write output to xml file.
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_nucnet,
    "%.15e"
  );

  Libnucnet__writeToXmlFile( p_my_nucnet, argv[4] );

  //============================================================================
  // Clean up.  Done.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
