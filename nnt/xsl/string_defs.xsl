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
  * \brief An xsl stylesheet to convert input into a header file.
  */

-->

<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="text"/>


<xsl:template match="/">

  <xsl:text>//! \file</xsl:text>
<xsl:text>
</xsl:text>
  <xsl:text>//! \brief An automatically generated file defining strings. To change it, edit xml/master.xml.</xsl:text>
<xsl:text>
</xsl:text>
<xsl:text>
</xsl:text>
<xsl:text>
</xsl:text>
  <xsl:text>#ifndef NNT_STRING_DEFS_H</xsl:text>
<xsl:text>
</xsl:text>
  <xsl:text>#define NNT_STRING_DEFS_H</xsl:text>
<xsl:text>
</xsl:text>
<xsl:text>
</xsl:text>
  <xsl:text>namespace nnt</xsl:text>
<xsl:text>
</xsl:text>
  <xsl:text>{</xsl:text>

<xsl:text>
</xsl:text>
<xsl:text>
</xsl:text>
   <xsl:for-each select="//string | //function">
     <xsl:sort select="key"/>
     <xsl:text>   const char </xsl:text>
     <xsl:value-of select="key"/>
     <xsl:text>[] = "</xsl:text>
     <xsl:value-of select="key_string"/>
     <xsl:text>";</xsl:text>
<xsl:text>
</xsl:text>
   </xsl:for-each>

<xsl:text>
</xsl:text>
<xsl:text>} // namespace nnt</xsl:text>
<xsl:text>
</xsl:text>
<xsl:text>
</xsl:text>
<xsl:text>#endif /* NNT_STRING_DEFS_H */</xsl:text>
<xsl:text>
</xsl:text>

</xsl:template>

</xsl:stylesheet>
