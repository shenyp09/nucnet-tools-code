///////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code to write out the abundance graph of a zone in dot
//!    format.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Include.
//##############################################################################

#include <fstream>
#include <ostream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>

#include <boost/graph/graphviz.hpp>

#include "nnt/iter.h"

#include "graph_helper.h"

int Counter;
std::ofstream parameter_out;

//##############################################################################
// create_graph().
//##############################################################################

void
create_data(
    nnt::Zone& zone,
    Libnucnet__NucView * p_nuc_view,
    double d_cutoff
)
{
    // std::string base = "dat_files";


    char s_output[100];
    snprintf(s_output, sizeof s_output, "dot_abund/abund_%05d.dat", Counter);

    std::ofstream my_out;
    my_out.open( s_output );
    Counter++;

    nnt::species_list_t species_list =
        nnt::make_species_list(
            Libnucnet__NucView__getNuc( p_nuc_view )
        );

    double a_abund_tmp[100];
    memset(a_abund_tmp, 0, 100 * sizeof(double));
    int max_a = 1;
    double d_xm = Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 );

    BOOST_FOREACH( nnt::Species species, species_list )
    {
        int i_a = Libnucnet__Species__getA( species.getNucnetSpecies() );
        if (i_a > max_a)
        {
            max_a = i_a;
        }

        double d_y =
            Libnucnet__Zone__getSpeciesAbundance(
                zone.getNucnetZone(),
                species.getNucnetSpecies()
            );
        a_abund_tmp[i_a] += d_y;
    }
    printf("%f\n", d_xm);
    for (int i = 1; i <= max_a; ++i)
    {
        char s[1000];
        snprintf(s, sizeof s, "%5d %10.2e %10.2e\n", i, a_abund_tmp[i] / d_xm, a_abund_tmp[i]*i);
        my_out << s;
    }

    my_out.close();

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

    Libnucnet *p_my_nucnet;
    Libnucnet__NucView * p_nuc_view;



    if ( argc != 6 )
    {
        fprintf(
            stderr,
            "\nUsage: %s xml_file nuc_xpath zone_xpath cutoff dot_file_base\n\n",
            argv[0]
        );
        fprintf(
            stderr,
            "  xml_file = network xml file\n\n"
        );
        fprintf(
            stderr,
            "  nuc_xpath = XPath expression to select nuclides\n\n"
        );
        fprintf(
            stderr,
            "  zone_xpath = XPath expression to select zones\n\n"
        );
        fprintf(
            stderr,
            "  cutoff = minimum abundance to include\n\n"
        );
        fprintf(
            stderr,
            "  dot_file_base = base name for output dot file\n\n"
        );
        return EXIT_FAILURE;
    }

    p_my_nucnet =
        Libnucnet__new_from_xml(
            argv[1],
            NULL,
            NULL,
            argv[3]
        );

    //============================================================================
    // Get subgraph view.
    //============================================================================

    p_nuc_view =
        Libnucnet__NucView__new(
            Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
            argv[2]
        );

    //============================================================================
    // Get zone list.
    //============================================================================

    Libnucnet__setZoneCompareFunction(
        p_my_nucnet,
        (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
    );

    nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

    //============================================================================
    // Iterate zones.
    //============================================================================

    BOOST_FOREACH( nnt::Zone zone, zone_list )
    {

        create_data( zone, p_nuc_view, atof( argv[4] ) );

    }

    //============================================================================
    // Clean up.
    //============================================================================

    Libnucnet__NucView__free( p_nuc_view );
    Libnucnet__free( p_my_nucnet );

    return EXIT_SUCCESS;

}
