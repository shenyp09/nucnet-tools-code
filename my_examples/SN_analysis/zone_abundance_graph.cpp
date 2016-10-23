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

#include <boost/graph/graphviz.hpp>

#include "nnt/iter.h"

#include "graph_helper.h"

#define D_BASE  120.
#define D_SCALE 10

Libnucnet__Nuc * p_nuc;

//##############################################################################
// The vertex writer.
//##############################################################################

struct
    my_vertex_writer
{
    my_vertex_writer( my_graph_t& g_ ) :           g( g_ ) {};
    template <class Vertex>
    void operator()( std::ostream& out, Vertex v )
    {
        Libnucnet__Species * p_sep, * p_species2;
        int i_z = Libnucnet__Species__getZ( g[v].getNucnetSpecies() );

        int i_n = Libnucnet__Species__getA( g[v].getNucnetSpecies() ) - i_z;

        int i_a = Libnucnet__Species__getA( g[v].getNucnetSpecies() );

        int i_z2 = i_z - 1;

        int i_a2 = i_a - 1;

        if ( ( i_z2 > 0 && i_a2 > 0 ) or ( i_z2 == 0 && i_a2 == 1 ) )
        {
            p_species2 =
                Libnucnet__Nuc__getSpeciesByZA(
                    p_nuc,
                    i_z2,
                    i_a2,
                    ""
                );
        }
        else
        {
            p_species2 = NULL;
        }

        double sep = 0.001;


        if ( p_species2 )
        {

            sep = 7.28897 +
                  Libnucnet__Species__getMassExcess( p_species2 ) -
                  Libnucnet__Species__getMassExcess( g[v].getNucnetSpecies() );

        }


        i_z2 = i_z;

        if ( ( i_z2 > 0 && i_a2 > 0 ) or ( i_z2 == 0 && i_a2 == 1 ) )
        {
            p_species2 =
                Libnucnet__Nuc__getSpeciesByZA(
                    p_nuc,
                    i_z2,
                    i_a2,
                    ""
                );
        }
        else
        {
            p_species2 = NULL;
        }

        double sen = 0.001;


        if ( p_species2 )
        {

            sen = 8.07132 +
                  Libnucnet__Species__getMassExcess( p_species2 ) -
                  Libnucnet__Species__getMassExcess( g[v].getNucnetSpecies() );

        }



        char * s_latex_string =
            Libnucnet__Species__createLatexString( g[v].getNucnetSpecies() );
        std::string style = "filled";
        double mex = Libnucnet__Species__getMassExcess(g[v].getNucnetSpecies());
        if (sep < 0 || sen < 0)
        {
            // printf("%s %f\n", Libnucnet__Species__getName(g[v].getNucnetSpecies()), mex);
            style = "invis";
        }
        out <<
            boost::format(
                " [texlbl=\"\" \
          label=\"\", \
          pos=\"%d,%d!\", \
          style=%s, fillcolor=\"%s\" \
         ]"
            ) %
            ( D_SCALE * i_n ) %
            ( D_SCALE * i_z ) %
            style.c_str() %
            g[v].getExtraData().color
            << std::endl;
        free( s_latex_string );
    }
    my_graph_t& g;
};

//##############################################################################
// The graph writer.
//##############################################################################

struct
    my_graph_writer
{
    my_graph_writer( nnt::Zone& zone_ ) : zone( zone_ ) {};
    void operator()( std::ostream& out ) const
    {
        out.precision(4);
        out << "graph [bgcolor=lightgrey  fontname=\"Helvetica\"];" << std::endl;
        out <<
            boost::format( "node [shape=box, color=\"black\", width=0.14, height=0.14, penwidth=1.];" )
            << std::endl;
        out <<
            boost::format(
                "label = \"time(s) = %g  T_9 = %g  rho(g/cc) = %g\" fontsize=\"30pt\" lp=\"400,800!\";"
            ) %
            zone.getProperty<double>( nnt::s_TIME ) %
            zone.getProperty<double>( nnt::s_T9 ) %
            zone.getProperty<double>( nnt::s_RHO ) <<
            std::endl;
    };
    nnt::Zone& zone;
};

//##############################################################################
// my_output().
//##############################################################################

void
my_output(
    my_graph_t& g,
    nnt::Zone zone,
    const char * s_output
)
{

    std::ofstream my_out;
    my_out.open( s_output );
    boost::write_graphviz(
        my_out,
        g,
        my_vertex_writer( g ),
        boost::default_writer(),
        my_graph_writer( zone )
    );
    my_out.close();

}

//##############################################################################
// create_graph().
//##############################################################################

my_graph_t
create_graph(
    nnt::Zone& zone,
    Libnucnet__NucView * p_nuc_view,
    double d_cutoff
)
{

    std::stringstream s_base, s_red, s_green, s_blue;
    std::string s_color, s_white = "#FFFFFF";
    my_graph_t g;

    my_graph_t::vertex_descriptor v;

    nnt::species_list_t species_list =
        nnt::make_species_list(
            Libnucnet__NucView__getNuc( p_nuc_view )
        );

    BOOST_FOREACH( nnt::Species species, species_list )
    {

        v = boost::add_vertex( g );

        g[v].setNucnetSpecies( species.getNucnetSpecies() );


        double d_y =
            Libnucnet__Zone__getSpeciesAbundance(
                zone.getNucnetZone(),
                species.getNucnetSpecies()
            );

        if ( d_y <= d_cutoff )
        {
            s_color = s_white;
        }
        else
        {
            s_red.str("");
            s_green.str("");
            s_blue.str("");
            s_red << std::hex <<
                  static_cast<int>(
                      254
                  );
            s_green << std::hex <<
                    static_cast<int>(
                        254 * (1 - log10(d_y) / log10(d_cutoff))
                    );
            s_blue << std::hex <<
                   static_cast<int>(
                       254 * (1 - log10(d_y) / log10(d_cutoff))
                   );

            s_color = "#" + s_red.str() + s_green.str() + s_blue.str();


            int red = 0, green = 0, blue = 0;
            double x = (log10(d_y) / log10(d_cutoff));
            red = (int)(488. - 233 * 4 * x);
            if (x <= 0.5)
            {
                green = (int)(22. + 223 * 4 * x);
            } else {
                green = (int)(22. + 233 * 4 * (1 - x));
            }
            blue = (int)(233 * 4 * x - 444.);
            if (red > 255)
            {
                red = 255;
            } else if (red < 22)
            {
                red = 22;
            }
            if (green > 255)
            {
                green = 255;
            } else if (green < 22)
            {
                green = 22;
            }
            if (blue > 255)
            {
                blue = 255;
            } else if (blue < 22)
            {
                blue = 22;
            }

            char hexcol[16];
            snprintf(hexcol, sizeof hexcol, "#%02x%02x%02x", red, green, blue);

            s_color = hexcol;
            // std::cout << log10(d_y) / log10(d_cutoff) << " " << s_color << std::endl;
        }

        g[v].updateExtraData( my_vertex_data( s_color, "box" ) );

    }

    return g;

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

    Libnucnet *p_my_nucnet;
    Libnucnet__NucView * p_nuc_view;


    p_nuc = Libnucnet__Nuc__new_from_xml( "../../data_pub/my_net.xml", NULL );

    my_graph_t g;

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

        g = create_graph( zone, p_nuc_view, atof( argv[4] ) );

        //==========================================================================
        // Add solar coloring.
        //==========================================================================

        // add_solar( g, p_nuc_view );

        //==========================================================================
        // Output.
        //==========================================================================

        std::stringstream ss;

        ss <<
           std::setw(5) <<
           std::setfill('0') <<
           Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 );

        std::string s_label(ss.str());

        std::string s_output =
            std::string( argv[5] ) +
            std::string( "_" ) +
            s_label;

        my_output( g, zone, s_output.c_str() );

    }

    //============================================================================
    // Clean up.
    //============================================================================

    Libnucnet__NucView__free( p_nuc_view );
    Libnucnet__free( p_my_nucnet );

    return EXIT_SUCCESS;

}
