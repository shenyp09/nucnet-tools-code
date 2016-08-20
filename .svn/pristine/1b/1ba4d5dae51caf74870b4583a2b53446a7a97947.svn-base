<?xml version="1.0"?>

<!--
   Copyright (c) 2015 Clemson University.
  
   This file was originally written by Bradley S. Meyer.
  
   This is free software; you can redistribute it and or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
  
   This software is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   USA
                                                                                
  /**
  * \file
  * \brief An xsl stylesheet to convert xml input into html documentation.
  */

  Stylesheet to display currently defined global strings and functions.  To
  process from nnt:

    xsltproc (-)(-)xinclude xsl/display.xsl xml/master.xml > display.html

  or, to select only subset of display,

    xsltproc (-)(-)xinclude (-)(-)param display_nodes NODES xsl/display.xsl xml/master.xml > display.html

  where NODES=//strings or NODES=//functions.

  Open display.html in a web browser.

-->

<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xi="http://www.w3.org/2001/XInclude"
  exclude-result-prefixes='xsl xi'>
<xsl:output method="xml" indent="yes"/>

<xsl:param name="display_nodes" select="//strings | //functions"/>

<xsl:template match="/">
  <xsl:apply-templates select = "$display_nodes" />
</xsl:template>

<xsl:template match="xi:include[@href][@parse='xml' or not(@parse)]">
 <xsl:apply-templates select="document(@href)" />
</xsl:template>

<xsl:template match="xi:include" />

<xsl:template match="functions">
  <html>
    <body>
      <h2>Functions</h2>
      <table border="1">
        <xsl:for-each select="function">
        <xsl:sort select="key"/>
          <tr>
            <td><b>Defined key: </b><xsl:value-of select="key"/><br/><b>Key string: </b><xsl:value-of select="key_string"/><br/><b>Description: </b><xsl:value-of select="doc"/><br/><b>Prototype: </b><xsl:value-of select="prototype"/></td>
          </tr>
        </xsl:for-each>
      </table>
    </body>
  </html>
</xsl:template>

<xsl:template match="strings">
  <html>
    <body>
      <h2>Strings</h2>
      <table border="1">
        <xsl:for-each select="string">
        <xsl:sort select="key"/>
          <tr>
            <td><b>Defined key: </b><xsl:value-of select="key"/><br/><b>Key string: </b><xsl:value-of select="key_string"/><br/><b>Description: </b><xsl:value-of select="doc"/></td>
          </tr>
        </xsl:for-each>
      </table>
    </body>
  </html>
</xsl:template>

</xsl:stylesheet>
