/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      Copyright (c) 2008-2016 Clemson University.
//
//      This distribution contains the source code
//      for the Clemson Webnucleo group's
//      libstatmech module, originally developed by Tianhong Yu and 
//      Bradley S. Meyer.
//
//      This is free software; you can redistribute it and/or modify it
//      under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This software is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this software (please see the "gnu_gpl.txt" file in the doc/
//      directory of this distribution); if not, write to the Free Software
//      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
//      USA
//
//      All libstatmech documentation is free documentation; permission is
//      granted to copy, distribute and/or modify the documentation under the
//      terms of the GNU Free Documentation License, Version 1.2 or any later
//      version published by the Free Software Foundation; with the Invariant
//      Sections, the Front-Cover Texts, and the Back-Cover Texts identified
//      in the documentation itself.  A copy of the license is included in the
//      file "gnu_fdl_v1.2.txt" in the doc/ directory in this distribution.
//   </license>
//
//   <description>
//     <abstract>
//       README file for the libstatmech src/ directory.
//     </abstract>
//     <keywords>
//       README, libstatmech, code, source
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="tyu" start_date="2008/05/20" />
//       <author userid="mbradle" start_date="2008/05/20" />
//     </current>
//     <previous>
//     </previous>
//   </authors>
//
// </file>
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////*/

################################################################################
# Overview of directory contents.
################################################################################

./README.txt

  Overview of directory contents and compatibility notes.

Libstatmech.c

  The libstatmech C code in the distribution for computing the statistical
  thermodynamics of fermion and boson gases.

Libstatmech.h

  The header file for Libstatmech.c.

################################################################################
# Dependencies.
################################################################################

libstatmech codes require libxml2 version 2.6.18 or later.  Earlier versions
of libxml2 may also work, but we haven't tested them.  To check whether
libxml2 is installed on a particular unix or linux system, type:

xml2-config --version

To install libxml2, please see http://www.libxml.org.

libstatmech codes require gsl version 1.9 or later.  Earlier versions
of gsl may also work, but we haven't tested them.  To check whether
gsl is installed on a particular unix or linux system, type:

gsl-config --version

To install gsl, please see http://www.gnu.org/software/gsl/

################################################################################
# Compatibility notes for the general user.
################################################################################

We've tested with gnu gcc 3.4.4 and Cygwin 1.5.24-2.

