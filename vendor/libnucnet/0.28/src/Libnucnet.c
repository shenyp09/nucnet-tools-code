/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      This file is part of the Clemson Webnucleo group's
//      libnucnet module, originally developed by David C. Adams
//      and Bradley S. Meyer.  For more information,
//      please see http://www.webnucleo.org.
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
//   </license>
//   <description>
//     <abstract>
//        Source code for the Libnucnet part of the libnucnet module.
//        For documentation, see the header file Libnucnet.h.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Standard library #include's.
//############################################################################*/

#include "Libnucnet.h"

/*##############################################################################
// Libnucnet__new().
//############################################################################*/

Libnucnet *Libnucnet__new( ) {

  Libnucnet *self;

  /*============================================================================
  // Get new network structure.
  //==========================================================================*/

  if( !( self = ( Libnucnet * ) malloc( sizeof( Libnucnet ) ) ) ) {

    LIBNUCNET__ERROR( "Virtual memory exhausted" );

  }

  /*============================================================================
  // Initialize structures.
  //==========================================================================*/

  self->pNet = Libnucnet__Net__new();

  /*============================================================================
  // Initialize Libnucnet quantities.
  //==========================================================================*/

  Libnucnet__initializeHash( self );
  self->pfZoneCompare = NULL;
  self->sxMassFractionFormat =
    xmlCharStrdup( ZONE_DEFAULT_MASS_FRACTION_FORMAT );

  return self;

}

/*##############################################################################
// Libnucnet__initializeHash()
//############################################################################*/

void Libnucnet__initializeHash( Libnucnet *self ) {

  self->pZoneHash = xmlHashCreate( 0 );

  return;

}
 

/*##############################################################################
// Libnucnet__new_from_xml().
//############################################################################*/

Libnucnet *
Libnucnet__new_from_xml(
  const char *s_xml_filename,
  const char *s_xpath_nuc,
  const char *s_xpath_reac,
  const char *s_xpath_zone
) {

  Libnucnet *self;

  /*============================================================================
  // Get new network structure.
  //==========================================================================*/

  if(
      !(
         self = Libnucnet__new()
      )
  ) {
    return NULL;
  }

  /*============================================================================
  // Read in and store nuclide and reaction data.
  //==========================================================================*/

  Libnucnet__Nuc__updateFromXml(
    self->pNet->pNuc,
    s_xml_filename,
    s_xpath_nuc
  );

  Libnucnet__Reac__updateFromXml(
    self->pNet->pReac,
    s_xml_filename,
    s_xpath_reac
  );

  /*============================================================================
  // Read in and store abundance data.
  //==========================================================================*/

  Libnucnet__assignZoneDataFromXml( self, s_xml_filename, s_xpath_zone );

  /*============================================================================
  // Done.  Return.
  //==========================================================================*/

  return self;

}

/*##############################################################################
// Libnucnet__assignZoneDataFromXml()
//############################################################################*/

void
Libnucnet__assignZoneDataFromXml(
  Libnucnet *self,
  const char *s_xml_filename,
  const char *s_xpath_suffix
) {

  xmlDocPtr p_doc;

  /*============================================================================
  // Create document
  //==========================================================================*/

  p_doc = xmlParseFile( s_xml_filename );

  if ( p_doc == NULL )
    LIBNUCNET__ERROR( "Document could not be opened" );

  if( xmlXIncludeProcess( p_doc ) == -1 )
    LIBNUCNET__ERROR( "Problem including xml" );

  /*============================================================================
  // Assign from XML document.
  //==========================================================================*/

  Libnucnet__assignZoneDataFromXmlDocument(
    self,
    p_doc,
    s_xpath_suffix
  );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// Libnucnet__assignZoneDataFromXmlDocument()
//############################################################################*/

void
Libnucnet__assignZoneDataFromXmlDocument(
  Libnucnet *self,
  xmlDocPtr p_doc,
  const char *s_xpath_suffix
) {

  xmlXPathContextPtr p_xpathCtx;
  xmlXPathObjectPtr p_xpathObj;
  xmlNodePtr p_zone;
  xmlChar *sx_xpath, *sx_xpath_suffix;
  xmlChar *sx_str, *sx_label_1, *sx_label_2, *sx_label_3;
  int i;

  Libnucnet__Zone *p_new_zone;

  /*============================================================================
  // Create context.
  //==========================================================================*/

  p_xpathCtx = xmlXPathNewContext( p_doc );

  /*============================================================================
  // Create main xpath expression.
  //==========================================================================*/

  sx_xpath = xmlCharStrdup( DOUBLE_SLASH );

  if( !sx_xpath ) {
     LIBNUCNET__ERROR( "Couldn't allocate memory for xpath string" );
  }

  sx_str = xmlCharStrdup( XPATH_ZONE_DATA );
  sx_xpath = xmlStrcat( sx_xpath, sx_str );
  xmlFree( sx_str );

  sx_xpath_suffix = xmlCharStrdup( s_xpath_suffix );

  if( sx_xpath_suffix ) {

     sx_xpath = xmlStrcat( sx_xpath, sx_xpath_suffix );

     if( !sx_xpath ) {
        LIBNUCNET__ERROR( "Couldn't allocate memory for xpath string" );
     }

     xmlFree( sx_xpath_suffix );

  }

  /*============================================================================
  // Loop over zones.
  //==========================================================================*/

  p_xpathObj = xmlXPathEvalExpression( sx_xpath, p_xpathCtx );
  xmlFree( sx_xpath );

  if( !p_xpathObj ) {
    LIBNUCNET__ERROR( "Couldn't evaluate xpath expression" );
  }

  if ( !p_xpathObj->nodesetval )
    LIBNUCNET__ERROR( "No zones in input" );

  for( i = 0; i < p_xpathObj->nodesetval->nodeNr; i++ ) {

    /*==========================================================================
    // Add zone.
    //========================================================================*/

    p_zone = p_xpathObj->nodesetval->nodeTab[i];

    sx_str = xmlCharStrdup( LABEL_1 );
    sx_label_1 = xmlGetProp( p_zone, sx_str );
    xmlFree( sx_str );

    sx_str = xmlCharStrdup( LABEL_2 );
    sx_label_2 = xmlGetProp( p_zone, sx_str );
    xmlFree( sx_str );

    sx_str = xmlCharStrdup( LABEL_3 );
    sx_label_3 = xmlGetProp( p_zone, sx_str );
    xmlFree( sx_str );

    p_new_zone =
      Libnucnet__Zone__new( 
        self->pNet,
        (const char *) sx_label_1,
        (const char *) sx_label_2,
        (const char *)sx_label_3
      );

    if( !p_new_zone ) LIBNUCNET__ERROR( "Couldn't create new zone" );

    /*==========================================================================
    // Free zone labels since new ones allocated in addZone.
    //========================================================================*/

    xmlFree( sx_label_1 );
    xmlFree( sx_label_2 );
    xmlFree( sx_label_3 );

    /*==========================================================================
    // Add properties and abundances.
    //========================================================================*/

    Libnucnet__Zone__assignOptionalPropertiesFromXml(
      p_new_zone,
      p_xpathObj->nodesetval->nodeTab[i]
    );

    Libnucnet__Zone__assignAbundancesFromXml(
      p_new_zone,
      p_xpathObj->nodesetval->nodeTab[i]
    );

    /*==========================================================================
    // Add zone.
    //========================================================================*/

    if( !Libnucnet__addZone( self, p_new_zone ) )
      LIBNUCNET__ERROR( "Couldn't add zone" );

  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlXPathFreeObject( p_xpathObj );
  xmlXPathFreeContext( p_xpathCtx );

}

/*##############################################################################
// Libnucnet__Zone__assignAbundancesFromXml()
//############################################################################*/

void
Libnucnet__Zone__assignAbundancesFromXml( 
  Libnucnet__Zone *self,
  xmlNodePtr p_zone_node
) {

  xmlXPathContextPtr p_xpathCtx;
  xmlXPathObjectPtr p_xpathObj;
  xmlNodePtr p_data;
  xmlChar *sx_str, *sx_z, *sx_a, *sx_m, *sx_name, *sx_prop;
  Libnucnet__Species *p_species;
  unsigned int i_z = 0, i_a = 0;
  int i;
  double d_x = 0.;

  /*============================================================================
  // Get context.
  //==========================================================================*/

  p_xpathCtx = xmlXPathNewContext( p_zone_node->doc );
  p_xpathCtx->node = p_zone_node;

  /*============================================================================
  // Loop over nuclide mass fractions.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XPATH_ZONE_MASS_FRACTIONS );
  p_xpathObj = xmlXPathEvalExpression( sx_str, p_xpathCtx );
  xmlFree( sx_str );

  if( !p_xpathObj ) {
    LIBNUCNET__ERROR( "Couldn' evaluate xpath expression" );
  }

  for( i = 0; i < p_xpathObj->nodesetval->nodeNr; i++ )
  {

    sx_name = xmlCharStrdup( SPECIES_NAME );

    sx_prop =
       xmlGetProp(
         p_xpathObj->nodesetval->nodeTab[i],
         sx_name
       );

    xmlFree( sx_name );

    p_data = p_xpathObj->nodesetval->nodeTab[i]->children;

    sx_z = xmlCharStrdup( ATOMIC_NUMBER );
    sx_a = xmlCharStrdup( MASS_NUMBER );
    sx_m = xmlCharStrdup( MASS_FRACTION );

    while( p_data ) {

      if( xmlStrEqual( sx_z, p_data->name ) )
          i_z =
            (unsigned int)
            xmlXPathCastStringToNumber( p_data->children->content );
      else if( xmlStrEqual( sx_a, p_data->name ) )
          i_a =
            (unsigned int)
            xmlXPathCastStringToNumber( p_data->children->content );
      else if( xmlStrEqual( sx_m, p_data->name ) )
          d_x =
            xmlXPathCastStringToNumber( p_data->children->content );

      p_data = p_data->next;

    }

    xmlFree( sx_z );
    xmlFree( sx_a );
    xmlFree( sx_m );

    if( sx_prop )
    {
      p_species = 
        Libnucnet__Nuc__getSpeciesByName(
          self->pNet->pNuc,
          (char *) sx_prop
        );
      if( !p_species )
      {
        fprintf( stderr, "Species: %s\n", (char *) sx_prop );
        LIBNUCNET__ERROR(
          "Can't assign initial mass fraction--species not present"
        );
      }
      if(
          Libnucnet__Species__getZ( p_species ) != i_z ||
          Libnucnet__Species__getA( p_species ) != i_a
      )
      {
        fprintf(
          stderr,
          "Zone: %s %s %s  Species: %s\n",
          Libnucnet__Zone__getLabel( self, 1 ),
          Libnucnet__Zone__getLabel( self, 2 ),
          Libnucnet__Zone__getLabel( self, 3 ),
          (char *) sx_prop
        );
        LIBNUCNET__ERROR(
          "Mismatch between species name and z or a"
        );
      }
      xmlFree( sx_prop );
    }
    else
    {
      p_species = 
        Libnucnet__Nuc__getSpeciesByZA( self->pNet->pNuc, i_z, i_a, "" );
      if ( !p_species ) {
        printf( "Species: Z = %u, A = %u\n", i_z, i_a );
        LIBNUCNET__ERROR(
          "Can't assign initial mass fraction--species not present"
        );
      } 
    }
    
    Libnucnet__Zone__updateSpeciesAbundance(
      self, p_species, d_x / p_species->iA
    );

  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlXPathFreeObject( p_xpathObj );
  xmlXPathFreeContext( p_xpathCtx );

}

/*##############################################################################
// Libnucnet__Zone__assignOptionalPropertiesFromXml()
//############################################################################*/

void
Libnucnet__Zone__assignOptionalPropertiesFromXml( 
  Libnucnet__Zone *self, xmlNodePtr p_zone_node
) {

  xmlXPathContextPtr p_xpathCtx;
  xmlXPathObjectPtr p_xpathObj;
  xmlChar *sx_str, *sx_prop_name, *sx_prop_tag1, *sx_prop_tag2, *sx_data;
  int i;

  /*============================================================================
  // Get context.
  //==========================================================================*/

  p_xpathCtx = xmlXPathNewContext( p_zone_node->doc );
  p_xpathCtx->node = p_zone_node;

  /*============================================================================
  // Loop over optional properties.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XPATH_OPTIONAL_PROPERTY );
  p_xpathObj = xmlXPathEvalExpression( sx_str, p_xpathCtx );
  xmlFree( sx_str );

  if( !p_xpathObj ) {
    LIBNUCNET__ERROR( "Couldn' evaluate xpath expression" );
  }

  for( i = 0; i < p_xpathObj->nodesetval->nodeNr; i++ )
  {

    sx_str = xmlCharStrdup( PROPERTY_NAME );

    sx_prop_name =
       xmlGetProp(
         p_xpathObj->nodesetval->nodeTab[i],
         sx_str
       );

    xmlFree( sx_str );

    sx_str = xmlCharStrdup( PROPERTY_TAG1 );

    sx_prop_tag1 =
       xmlGetProp(
         p_xpathObj->nodesetval->nodeTab[i],
         sx_str
       );

    xmlFree( sx_str );

    sx_str = xmlCharStrdup( PROPERTY_TAG2 );

    sx_prop_tag2 =
       xmlGetProp(
         p_xpathObj->nodesetval->nodeTab[i],
         sx_str
       );

    xmlFree( sx_str );

    sx_data = xmlNodeGetContent( p_xpathObj->nodesetval->nodeTab[i] );

    xmlHashUpdateEntry3(
      self->pPropertyHash,
      sx_prop_name,
      sx_prop_tag1,
      sx_prop_tag2,
      sx_data,
      (xmlHashDeallocator) xmlFree
    ); 

    xmlFree( sx_prop_name );
    xmlFree( sx_prop_tag1 );
    xmlFree( sx_prop_tag2 );

  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlXPathFreeObject( p_xpathObj );
  xmlXPathFreeContext( p_xpathCtx );

}

/*##############################################################################
// Libnucnet__Zone__new()
//############################################################################*/

Libnucnet__Zone *
Libnucnet__Zone__new(
  Libnucnet__Net *p_net,
  const char *s_label1,
  const char *s_label2,
  const char *s_label3
)
{

  Libnucnet__Zone *p_zone;

  p_zone = ( Libnucnet__Zone * ) malloc( sizeof( Libnucnet__Zone ) );

  /*============================================================================
  // Link zone to network.
  //==========================================================================*/

  p_zone->pNet = p_net;

  /*============================================================================
  // Create hashes.
  //==========================================================================*/

  p_zone->pAbundanceHash = xmlHashCreate( 0 );
  p_zone->pAbundanceChangeHash = xmlHashCreate( 0 );
  p_zone->pRatesHash = xmlHashCreate( 0 );
  p_zone->pPropertyHash = xmlHashCreate( 0 );
  p_zone->pUserDataHash = xmlHashCreate( 0 );
  p_zone->pViewHash = xmlHashCreate( 0 );

  /*============================================================================
  // Initialize data.
  //==========================================================================*/

  p_zone->iComputeReverseDetailedBalance = 1;

  p_zone->pfScreeningFunction = NULL;
  p_zone->pScreeningData = NULL;
  p_zone->pfNseCorrectionFactorFunction = NULL;
  p_zone->pNseCorrectionFactorData = NULL;

  /*============================================================================
  // Assign zone labels.
  //==========================================================================*/

  if( s_label1 )
    p_zone->sxLabel[0] =
      xmlCharStrdup( s_label1 );
  else
    p_zone->sxLabel[0] =
      xmlCharStrdup( ZERO );

  if( s_label2 )
    p_zone->sxLabel[1] =
      xmlCharStrdup( s_label2 );
  else
    p_zone->sxLabel[1] =
      xmlCharStrdup( ZERO );

  if( s_label3 )
    p_zone->sxLabel[2] =
      xmlCharStrdup( s_label3 );
  else
    p_zone->sxLabel[2] =
      xmlCharStrdup( ZERO );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return p_zone;

}

/*##############################################################################
// Libnucnet__addZone()
//############################################################################*/

int
Libnucnet__addZone(
  Libnucnet *self, 
  Libnucnet__Zone *p_zone
)
{

  /*============================================================================
  // Add zone to hash.
  //==========================================================================*/

  if( !self->pZoneHash ) Libnucnet__initializeHash( self );

  if(
    xmlHashAddEntry3(
      self->pZoneHash,
      p_zone->sxLabel[0],
      p_zone->sxLabel[1],
      p_zone->sxLabel[2],
      p_zone
    )
  )
  {

    Libnucnet__Zone__hashFree( p_zone, p_zone->sxLabel[0] );
    return 0;

  }

  return 1;

}

/*##############################################################################
// Libnucnet__Zone__setScreeningFunction()
//############################################################################*/

void
Libnucnet__Zone__setScreeningFunction(
  Libnucnet__Zone *self,
  Libnucnet__Zone__screeningFunction pfScreeningFunction,
  void *p_user_data
)
{
  if( !self ) {
    LIBNUCNET__ERROR( "Invalid input" );
  }

  self->pfScreeningFunction = pfScreeningFunction;
  self->pScreeningData = p_user_data;

}

/*##############################################################################
// Libnucnet__Zone__getScreeningFunction().
//############################################################################*/

Libnucnet__Zone__screeningFunction
Libnucnet__Zone__getScreeningFunction(
  Libnucnet__Zone * self
)
{

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  return self->pfScreeningFunction;

}

/*##############################################################################
// Libnucnet__Zone__getScreeningData().
//############################################################################*/

void *
Libnucnet__Zone__getScreeningData(
  Libnucnet__Zone * self
)
{

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  return self->pScreeningData;

}

/*##############################################################################
// Libnucnet__Zone__setNseCorrectionFactorFunction()
//############################################################################*/

void
Libnucnet__Zone__setNseCorrectionFactorFunction(
  Libnucnet__Zone *self,
  Libnucnet__Species__nseCorrectionFactorFunction
    pfNseCorrectionFactorFunction,
  void *p_user_data
)
{
  if( !self ) {
    LIBNUCNET__ERROR( "Invalid input" );
  }

  self->pfNseCorrectionFactorFunction =
    pfNseCorrectionFactorFunction;
  self->pNseCorrectionFactorData = p_user_data;

}

/*##############################################################################
// Libnucnet__Zone__getNseCorrectionFactorFunction().
//############################################################################*/

Libnucnet__Species__nseCorrectionFactorFunction
Libnucnet__Zone__getNseCorrectionFactorFunction(
  Libnucnet__Zone * self
)
{

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  return self->pfNseCorrectionFactorFunction;

}

/*##############################################################################
// Libnucnet__Zone__getNseCorrectionFactorData().
//############################################################################*/

void *
Libnucnet__Zone__getNseCorrectionFactorData(
  Libnucnet__Zone * self
)
{

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  return self->pNseCorrectionFactorData;

}

/*##############################################################################
// Libnucnet__Net__new().
//############################################################################*/

Libnucnet__Net *Libnucnet__Net__new( ) {

  /*============================================================================
  // Get new network structure.
  //==========================================================================*/

  Libnucnet__Net *self = (Libnucnet__Net *) malloc( sizeof( Libnucnet__Net ) );

  /*============================================================================
  // Initialize structures.
  //==========================================================================*/

  self->pNuc = Libnucnet__Nuc__new();
  self->pReac = Libnucnet__Reac__new();

  /*============================================================================
  // Return.
  //==========================================================================*/

  return self;

}

/*##############################################################################
// Libnucnet__Net__new_from_xml().
//############################################################################*/

Libnucnet__Net *
Libnucnet__Net__new_from_xml(
  const char *s_xml_filename, const char *s_xpath_nuc, const char *s_xpath_reac
) {

  Libnucnet__Net *self;

  /*============================================================================
  // Allocate memory for new network structure.
  //==========================================================================*/

  self = ( Libnucnet__Net * ) malloc( sizeof( Libnucnet__Net ) );

  if( !self ) {
    return NULL;
  }

  /*============================================================================
  // Get nuclear data.
  //==========================================================================*/

  if(
       !(
          self->pNuc =
            Libnucnet__Nuc__new_from_xml( s_xml_filename, s_xpath_nuc )
       )
    )
  {

    return NULL;

  }

  /*============================================================================
  // Get reaction structure and store data.
  //==========================================================================*/

  if(  
       !(
          self->pReac =
            Libnucnet__Reac__new_from_xml( s_xml_filename, s_xpath_reac )
       ) 
    )
  {

    return NULL;

  }

  /*============================================================================
  // Return.
  //==========================================================================*/

  return self;

}

/*##############################################################################
// Libnucnet__Net__updateFromXml().
//############################################################################*/

void
Libnucnet__Net__updateFromXml(
  Libnucnet__Net *self,
  const char *s_xml_filename,
  const char *s_xpath_nuc,
  const char *s_xpath_reac
) {

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Update nuclear data.
  //==========================================================================*/

  Libnucnet__Nuc__updateFromXml( self->pNuc, s_xml_filename, s_xpath_nuc );

  /*============================================================================
  // Get reaction structure and store data.
  //==========================================================================*/

  Libnucnet__Reac__updateFromXml( self->pReac, s_xml_filename, s_xpath_reac );

  /*============================================================================
  // Return.
  //==========================================================================*/

}

/*##############################################################################
// Libnucnet__NetView__new().
//############################################################################*/

Libnucnet__NetView *
Libnucnet__NetView__new(
  Libnucnet__Net * self,
  const char * s_nuc_xpath,
  const char * s_reac_xpath
)
{

  Libnucnet__NetView * p_view;
  Libnucnet__NucView * p_nuc_view;
  Libnucnet__ReacView * p_reac_view;

  p_view =
    ( Libnucnet__NetView * ) malloc( sizeof( Libnucnet__NetView ) );

  if( !p_view ) LIBNUCNET__ERROR( "Couldn't allocate memory for view" );

  p_view->pNet =
    ( Libnucnet__Net * ) malloc( sizeof( Libnucnet__Net ) );

  p_view->pParentNet = self;

  p_nuc_view = Libnucnet__NucView__new( self->pNuc, s_nuc_xpath );

  p_view->pNet->pNuc = p_nuc_view->pNuc;

  free( p_nuc_view );

  p_view->pNet->pReac =
    Libnucnet__Reac__createNewReacForView( self->pReac );

  p_reac_view = Libnucnet__ReacView__new( self->pReac, s_reac_xpath );

  xmlHashScan(
    p_reac_view->pReac->pReactionHash,
    (xmlHashScanner) Libnucnet__NetView__set_valid_reactions,
    p_view
  );

  Libnucnet__ReacView__free( p_reac_view );

  return p_view;

}
  
/*##############################################################################
// Libnucnet__NetView__getNet().
//############################################################################*/

Libnucnet__Net *
Libnucnet__NetView__getNet(
  Libnucnet__NetView * self
)
{

  if( !self) LIBNUCNET__ERROR( "Invalid input" );

  return self->pNet;

}

/*##############################################################################
// Libnucnet__NetView__copy().
//############################################################################*/

Libnucnet__NetView *
Libnucnet__NetView__copy(
  Libnucnet__NetView * self
)
{

  Libnucnet__Nuc * p_nuc;
  Libnucnet__Reac * p_reac;
  Libnucnet__NetView * p_new_view;

  if( !self ) LIBNUCNET__ERROR( "Copying invalid net view" );

  p_new_view =
    ( Libnucnet__NetView * ) malloc( sizeof( Libnucnet__NetView ) );

  if( !p_new_view ) LIBNUCNET__ERROR( "Couldn't allocate memory for view" );

  p_new_view->pNet =
    ( Libnucnet__Net * ) malloc( sizeof( Libnucnet__Net ) );

  if( !p_new_view ) LIBNUCNET__ERROR( "Couldn't allocate memory for view" );

  p_new_view->pParentNet = self->pParentNet;

  p_nuc = 
    Libnucnet__Net__getNuc(
      Libnucnet__NetView__getNet( self )
    );

  p_new_view->pNet->pNuc =
    Libnucnet__Nuc__createNewNucForView( p_nuc );

  xmlHashScan(
    p_nuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__NetView__hash_element_copier,
    p_new_view->pNet->pNuc->pSpeciesHash
  );

  p_reac =
    Libnucnet__Net__getReac(
      Libnucnet__NetView__getNet( self )
    );

  p_new_view->pNet->pReac =
    Libnucnet__Reac__createNewReacForView( p_reac );

  xmlHashScan(
    p_reac->pReactionHash,
    (xmlHashScanner) Libnucnet__NetView__hash_element_copier,
    p_new_view->pNet->pReac->pReactionHash
  );

  return p_new_view;

}

/*##############################################################################
// Libnucnet__NetView__hash_element_copier().
//############################################################################*/

void
Libnucnet__NetView__hash_element_copier(
  void *p_data,
  xmlHashTablePtr p_hash,
  xmlChar *sx_name
)
{

  xmlHashAddEntry( p_hash, sx_name, p_data );

} 

/*##############################################################################
// Libnucnet__NetView__set_valid_reactions().
//############################################################################*/

void
Libnucnet__NetView__set_valid_reactions(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__NetView * p_view,
  xmlChar * sx_string
)
{

  if( !sx_string ) LIBNUCNET__ERROR( "Invalid reaction" );

  if(
      Libnucnet__Net__isValidReaction(
        Libnucnet__NetView__getNet( p_view ),
        p_reaction
      )
  ) 
    Libnucnet__NetView__addReaction(
      p_view,
      p_reaction
    );

} 

/*##############################################################################
// Libnucnet__Net__getNumberOfValidReactions()
//############################################################################*/

size_t
Libnucnet__Net__getNumberOfValidReactions( Libnucnet__Net *self )
{

  size_t i_result;
  Libnucnet__NetView * p_view;

  p_view = Libnucnet__NetView__new( self, NET_EMPTY_STRING, NET_EMPTY_STRING );

  i_result = (size_t) xmlHashSize( p_view->pNet->pReac->pReactionHash );

  Libnucnet__NetView__free( p_view );

  return i_result;

}

/*##############################################################################
// Libnucnet__Net__checkReactionForBaryonConservation().
//############################################################################*/

void
Libnucnet__Net__checkReactionForBaryonConservation(
  const Libnucnet__Net *self,
  const Libnucnet__Reaction *p_reaction
) {

  typedef struct {
    int iBaryonNumberChange;
    const Libnucnet__Net *pNet;
  } user_data;

  user_data *p_user_data;

  /*============================================================================
  // Allocate memory.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data )
    LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_user_data->pNet = self;
  p_user_data->iBaryonNumberChange = 0;

  /*============================================================================
  // Compute change from reactants.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pReactantList,
    (xmlListWalker) Libnucnet__Net__checkReactionForBaryonConservationWalker,
    p_user_data
  );

  p_user_data->iBaryonNumberChange *= -1;

  /*============================================================================
  // Compute change from products.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pProductList,
    (xmlListWalker) Libnucnet__Net__checkReactionForBaryonConservationWalker,
    p_user_data
  );

  /*============================================================================
  // Check if baryon number change is not zero.
  //==========================================================================*/

  if( p_user_data->iBaryonNumberChange ) {
    fprintf( stderr, "%s\n", Libnucnet__Reaction__getString( p_reaction ) );
    LIBNUCNET__ERROR( "Baryon number not conserved" );
  }

  free( p_user_data );

}

/*##############################################################################
// Libnucnet__Net__checkReactionForBaryonConservationWalker().
//############################################################################*/

int
Libnucnet__Net__checkReactionForBaryonConservationWalker(
  Libnucnet__Reaction__Element *p_element, void *p_data
)
{

  typedef struct {
    int iBaryonNumberChange;
    Libnucnet__Net *pNet;
  } user_data;

  user_data *p_user_data = ( user_data * ) p_data;

  p_user_data->iBaryonNumberChange +=
    (int)
    Libnucnet__Species__getA(
      Libnucnet__Nuc__getSpeciesByName(
        p_user_data->pNet->pNuc,
        (char *) p_element->sxName
      )
    );

  return 1;

}

/*##############################################################################
// Libnucnet__Net__checkReactionForChargeConservation().
//############################################################################*/

void
Libnucnet__Net__checkReactionForChargeConservation(
  const Libnucnet__Net *self,
  const Libnucnet__Reaction *p_reaction
) {

  typedef struct {
    const Libnucnet__Net *pNet;
    int iChargeChange;
  } user_data;

  user_data *p_user_data;

  /*============================================================================
  // Allocate memory.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_user_data->pNet = self;
  p_user_data->iChargeChange = 0;

  /*============================================================================
  // Compute change from reactants.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pReactantList,
    (xmlListWalker) Libnucnet__Net__checkReactionForChargeConservationWalker,
    p_user_data
  );

  xmlListWalk(
    p_reaction->pOtherReactantList,
    (xmlListWalker) Libnucnet__Net__checkReactionForChargeConservationWalker,
    p_user_data
  );

  p_user_data->iChargeChange *= -1;

  /*============================================================================
  // Compute change from products.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pProductList,
    (xmlListWalker) Libnucnet__Net__checkReactionForChargeConservationWalker,
    p_user_data
  );

  xmlListWalk(
    p_reaction->pOtherProductList,
    (xmlListWalker) Libnucnet__Net__checkReactionForChargeConservationWalker,
    p_user_data
  );
  
  /*============================================================================
  // Check if charge change is not zero.
  //==========================================================================*/

  if( p_user_data->iChargeChange ) {
    fprintf( stderr, "%s\n", Libnucnet__Reaction__getString( p_reaction ) );
    LIBNUCNET__ERROR( "Charge not conserved" );
  }

  free( p_user_data );

}

/*##############################################################################
// Libnucnet__Net__checkReactionForChargeConservationWalker().
//############################################################################*/

int
Libnucnet__Net__checkReactionForChargeConservationWalker(
  Libnucnet__Reaction__Element *p_element, void *p_data
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    int iChargeChange;
  } user_data;

  user_data *p_user_data = ( user_data * ) p_data;

  if( Libnucnet__Reaction__Element__isNuclide( p_element ) ) {

    p_user_data->iChargeChange +=
      (int)
      Libnucnet__Species__getZ( 
        Libnucnet__Nuc__getSpeciesByName(
          p_user_data->pNet->pNuc,
          (char *) p_element->sxName
        )
      );

  } else {

    if( strcmp( (char *) p_element->sxName, ELECTRON ) == 0 ) {
      p_user_data->iChargeChange--;
    }
    else if( strcmp( (char *) p_element->sxName, POSITRON ) == 0 ) {
      p_user_data->iChargeChange++;
    }
    else if( strcmp( (char *) p_element->sxName, MU ) == 0 ) {
      p_user_data->iChargeChange--;
    }
    else if( strcmp( (char *) p_element->sxName, ANTIMU ) == 0 ) {
      p_user_data->iChargeChange++;
    }
    else if( strcmp( (char *) p_element->sxName, TAU ) == 0 ) {
      p_user_data->iChargeChange--;
    }
    else if( strcmp( (char *) p_element->sxName, ANTITAU ) == 0 ) {
      p_user_data->iChargeChange++;
    }

  }

  return 1;

}

/*##############################################################################
// Libnucnet__Net__checkReactionForLeptonConservation().
//############################################################################*/

void
Libnucnet__Net__checkReactionForLeptonConservation(
  const Libnucnet__Net *self,
  const Libnucnet__Reaction *p_reaction
) {

  typedef struct {
    const Libnucnet__Net *pNet;
    int iElectronLeptonNumberChange;
    int iMuLeptonNumberChange;
    int iTauLeptonNumberChange;
  } user_data;

  user_data *p_user_data;

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid network structure" );
  }

  /*============================================================================
  // Allocate memory.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_user_data->pNet = self;
  p_user_data->iElectronLeptonNumberChange = 0;
  p_user_data->iMuLeptonNumberChange = 0;
  p_user_data->iTauLeptonNumberChange = 0;

  /*============================================================================
  // Check reactants.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pOtherReactantList,
    (xmlListWalker) Libnucnet__Net__checkReactionForLeptonConservationWalker,
    p_user_data
  );

  p_user_data->iElectronLeptonNumberChange *= -1;
  p_user_data->iMuLeptonNumberChange *= -1;
  p_user_data->iTauLeptonNumberChange *= -1;

  /*============================================================================
  // Check products.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pOtherProductList,
    (xmlListWalker) Libnucnet__Net__checkReactionForLeptonConservationWalker,
    p_user_data
  );

  /*============================================================================
  // Check if lepton numbers not conserved.
  //==========================================================================*/

  if( p_user_data->iElectronLeptonNumberChange ) {
    fprintf( stderr, "%s\n", Libnucnet__Reaction__getString( p_reaction ) );
    LIBNUCNET__ERROR( "Electron lepton number not conserved" );
  }

  if( p_user_data->iMuLeptonNumberChange ) {
    fprintf( stderr, "%s\n", Libnucnet__Reaction__getString( p_reaction ) );
    LIBNUCNET__ERROR( "Muon lepton number not conserved" );
  }

  if( p_user_data->iTauLeptonNumberChange ) {
    fprintf( stderr, "%s\n", Libnucnet__Reaction__getString( p_reaction ) );
    LIBNUCNET__ERROR( "Tauon lepton number not conserved" );
  }

  /*============================================================================
  // Done.  Free user data.
  //==========================================================================*/

  free( p_user_data );

}

/*##############################################################################
// Libnucnet__Net__checkReactionForLeptonConservationWalker().
//############################################################################*/

int
Libnucnet__Net__checkReactionForLeptonConservationWalker(
  Libnucnet__Reaction__Element *p_element, void *p_data
) {

  typedef struct {
    Libnucnet__Net *pNet;
    int iElectronLeptonNumberChange;
    int iMuLeptonNumberChange;
    int iTauLeptonNumberChange;
  } user_data;

  user_data *p_user_data = ( user_data * ) p_data;

  /*============================================================================
  // Check electron lepton number change.
  //==========================================================================*/

  if( strcmp( (char *) p_element->sxName, ELECTRON_NEUTRINO ) == 0 ) {
    p_user_data->iElectronLeptonNumberChange++;
  }
  else if(
    strcmp( (char *) p_element->sxName, ELECTRON_ANTINEUTRINO ) == 0
  ) {
    p_user_data->iElectronLeptonNumberChange--;
  }
  else if( strcmp( (char *) p_element->sxName, ELECTRON ) == 0 ) {
    p_user_data->iElectronLeptonNumberChange++;
  }
  else if( strcmp( (char *) p_element->sxName, POSITRON ) == 0 ) {
    p_user_data->iElectronLeptonNumberChange--;
  }

  /*============================================================================
  // Check muon lepton number change.
  //==========================================================================*/

  if( strcmp( (char *) p_element->sxName, MU_NEUTRINO ) == 0 ) {
    p_user_data->iMuLeptonNumberChange++;
  }
  else if(
    strcmp( (char *) p_element->sxName, MU_ANTINEUTRINO ) == 0
  ) {
    p_user_data->iMuLeptonNumberChange--;
  }
  else if( strcmp( (char *) p_element->sxName, MU ) == 0 ) {
    p_user_data->iMuLeptonNumberChange++;
  }
  else if( strcmp( (char *) p_element->sxName, ANTIMU ) == 0 ) {
    p_user_data->iMuLeptonNumberChange--;
  }

  /*============================================================================
  // Check tauon lepton number change.
  //==========================================================================*/

  if( strcmp( (char *) p_element->sxName, TAU_NEUTRINO ) == 0 ) {
    p_user_data->iTauLeptonNumberChange++;
  }
  else if(
    strcmp( (char *) p_element->sxName, TAU_ANTINEUTRINO ) == 0
  ) {
    p_user_data->iTauLeptonNumberChange--;
  }
  else if( strcmp( (char *) p_element->sxName, TAU ) == 0 ) {
    p_user_data->iTauLeptonNumberChange++;
  }
  else if( strcmp( (char *) p_element->sxName, ANTITAU ) == 0 ) {
    p_user_data->iTauLeptonNumberChange--;
  }

  return 1;

}

/*##############################################################################
// Libnucnet__Zone__computeRates().
//############################################################################*/

void
Libnucnet__Zone__computeRates(
  Libnucnet__Zone *self,
  double d_t9,
  double d_rho
)
{

  typedef struct {
    const Libnucnet__Zone *pZone;
    double dT9;
    double dRho;
    double dYe;
  } extra_data;

  extra_data *p_extra_data;
  Libnucnet__NetView * p_view;
    
  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !self ) LIBNUCNET__ERROR( "Invalid input" );

  /*============================================================================
  // Get valid view.
  //==========================================================================*/

  p_view = Libnucnet__Zone__getEvolutionNetView( self );

  /*============================================================================
  // Set up work structure.
  //==========================================================================*/

  p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  p_extra_data->dT9 = d_t9;
  p_extra_data->dRho = d_rho;
  p_extra_data->dYe = Libnucnet__Zone__computeZMoment( self, 1 );

  p_extra_data->pZone = self;

  /*============================================================================
  // Compute rates.
  //==========================================================================*/

  xmlHashScan(
    p_view->pNet->pReac->pReactionHash,
    (xmlHashScanner) Libnucnet__Zone__computeRatesCallback,
    p_extra_data
  );

  /*============================================================================
  // Apply screening.
  //==========================================================================*/

  if( self->pfScreeningFunction ) {

    xmlHashScan(
      p_view->pNet->pReac->pReactionHash,
      (xmlHashScanner) Libnucnet__Zone__applyScreeningCallback,
      p_extra_data
    );

  }

  /*============================================================================
  // Clean up.  Done.
  //==========================================================================*/

  free( p_extra_data );

}

/*##############################################################################
// Libnucnet__Zone__computeRatesCallback().
//############################################################################*/

void
Libnucnet__Zone__computeRatesCallback(
  Libnucnet__Reaction *p_reaction,
  void *p_data,
  xmlChar *sx_reaction
) {

  typedef struct {
    Libnucnet__Zone *pZone;
    double dT9;
    double dRho;
    double dYe;
  } extra_data;
    
  Libnucnet__Rates *p_rates;
  double d_forward, d_reverse;
  extra_data *p_extra_data;
  size_t i_number_nuclide_reactants, i_number_nuclide_products;

  /*==========================================================================
  // Check input.
  //========================================================================*/

  p_extra_data = ( extra_data * ) p_data;

  /*==========================================================================
  // Compute rates.
  //========================================================================*/

  if(
     Libnucnet__Zone__isComputingReverseRatesFromDetailedBalance(
       p_extra_data->pZone
     )
  )
  {

    Libnucnet__Net__computeRatesForReaction(
      p_extra_data->pZone->pNet,
      p_reaction,
      p_extra_data->dT9,
      p_extra_data->dRho,
      Libnucnet__Zone__getDataForUserRateFunction(
	p_extra_data->pZone,
	Libnucnet__Reaction__getRateFunctionKey( p_reaction )
      ),
      &d_forward,
      &d_reverse
    );

  }
  else
  {

    d_forward =
      Libnucnet__Reaction__computeRate(
        p_reaction,
        p_extra_data->dT9,
        Libnucnet__Zone__getDataForUserRateFunction(
          p_extra_data->pZone,
          Libnucnet__Reaction__getRateFunctionKey( p_reaction )
        )
      );

    d_reverse = 0;

  }
    
  /*==========================================================================
  // Create Rates structure.
  //========================================================================*/

  p_rates = (Libnucnet__Rates *) malloc( sizeof( Libnucnet__Rates ) );

  /*==========================================================================
  // Count reactants and products.
  //========================================================================*/

  i_number_nuclide_reactants =
    (size_t) xmlListSize( p_reaction->pReactantList );

  i_number_nuclide_products = (size_t) xmlListSize( p_reaction->pProductList );

  /*==========================================================================
  // Assign density term.
  //========================================================================*/

  p_rates->dForward =
    d_forward *
    pow(
      p_extra_data->dRho,
      ( double ) i_number_nuclide_reactants - 1.
    );

  p_rates->dReverse =
    d_reverse *
    pow(
      p_extra_data->dRho,
      ( double ) i_number_nuclide_products - 1.
    );

  /*==========================================================================
  // Assign rates to hash.
  //========================================================================*/

  if( 
       (
          xmlHashUpdateEntry(
            p_extra_data->pZone->pRatesHash,
            sx_reaction,
            p_rates,
            (xmlHashDeallocator) Libnucnet__Rates__free
          )
       )
  ) {
      LIBNUCNET__ERROR( "Couldn't update entry" );
  }

}

/*##############################################################################
// Libnucnet__Zone__applyScreeningCallback().
//############################################################################*/

void
Libnucnet__Zone__applyScreeningCallback(
  Libnucnet__Reaction *p_reaction,
  void *p_data,
  xmlChar *sx_reaction
)
{

  typedef struct {
    Libnucnet__Zone *pZone;
    double dT9;
    double dRho;
    double dYe;
  } extra_data;
    
  extra_data *p_extra_data;
  double d_forward, d_reverse;

  if( !sx_reaction ) LIBNUCNET__ERROR( "Invalid input" );

  /*==========================================================================
  // Get work structure.
  //========================================================================*/

  p_extra_data = ( extra_data * ) p_data;

  /*==========================================================================
  // Get the rates.
  //========================================================================*/

  Libnucnet__Zone__getRatesForReaction(
    p_extra_data->pZone,
    p_reaction,
    &d_forward,
    &d_reverse
  );

  /*==========================================================================
  // Apply screening.
  //========================================================================*/

  p_extra_data->pZone->pfScreeningFunction(
    p_extra_data->pZone,
    p_reaction,
    p_extra_data->dT9,
    p_extra_data->dRho,
    p_extra_data->dYe,
    &d_forward,
    &d_reverse
  );

  Libnucnet__Zone__updateRatesForReaction(
    p_extra_data->pZone,
    p_reaction,
    d_forward,
    d_reverse
  );

}

/*##############################################################################
// Libnucnet__Net__computeReverseRatioCorrectionFactorForReaction().
//############################################################################*/

double
Libnucnet__Net__computeReverseRatioCorrectionFactorForReaction(
  Libnucnet__Net *self,
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  double d_rho,
  double d_ye,
  Libnucnet__Species__nseCorrectionFactorFunction pf_correction_function,
  void *p_user_data
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    double dCorr;
    double dT9;
    double dRho;
    double dYe;
    Libnucnet__Species__nseCorrectionFactorFunction pfFunc;
    void *pUserData;
  } extra_data;

  double d_corr;
  extra_data *p_extra_data;

  if( !self || !p_reaction )
    LIBNUCNET__ERROR( "Invalid input" );

  /*==========================================================================
  // Allocate memory for extra data.
  //========================================================================*/

  p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_extra_data ) LIBNUCNET__ERROR( "Couldn't allocate memory" );

  /*==========================================================================
  // Assign data to structure.
  //========================================================================*/

  p_extra_data->pNet = self;
  p_extra_data->dCorr = 0.;
  p_extra_data->dT9 = d_t9;
  p_extra_data->dRho = d_rho;
  p_extra_data->dYe = d_ye;
  p_extra_data->pfFunc = pf_correction_function;
  p_extra_data->pUserData = p_user_data;

  /*==========================================================================
  // Iterate reactants and products.
  //========================================================================*/

  if( pf_correction_function )
  {

   /*----------------------------------------------------------------------
   // First reactants.
   //--------------------------------------------------------------------*/

   xmlListWalk(
     p_reaction->pReactantList,
     (xmlListWalker)
        Libnucnet__Net__computeReverseRatioCorrectionFactorForReactionWalker,
     p_extra_data
   );

   p_extra_data->dCorr *= -1.;

   /*----------------------------------------------------------------------
   // Now products.
   //--------------------------------------------------------------------*/

   xmlListWalk(
     p_reaction->pProductList,
     (xmlListWalker)
        Libnucnet__Net__computeReverseRatioCorrectionFactorForReactionWalker,
     p_extra_data
   );

  }

  d_corr = exp( p_extra_data->dCorr );

  free( p_extra_data );

  return d_corr;

}

/*##############################################################################
// Libnucnet__Net__computeReverseRatioCorrectionFactorForReactionWalker().
//############################################################################*/

int
Libnucnet__Net__computeReverseRatioCorrectionFactorForReactionWalker(
  Libnucnet__Reaction__Element *p_element, void *p_data
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    double dCorr;
    double dT9;
    double dRho;
    double dYe;
    Libnucnet__Species__nseCorrectionFactorFunction pfFunc;
    void *pUserData;
  } extra_data;

  extra_data *p_extra_data = ( extra_data * ) p_data;

  p_extra_data->dCorr -=
    p_extra_data->pfFunc(
      Libnucnet__Nuc__getSpeciesByName(
        p_extra_data->pNet->pNuc,
        (char *) p_element->sxName
      ),
      p_extra_data->dT9,
      p_extra_data->dRho,
      p_extra_data->dYe,
      p_extra_data->pUserData
    ); 

  return 1;

}

/*##############################################################################
// Libnucnet__Net__computeRatesForReaction().
//############################################################################*/

void
Libnucnet__Net__computeRatesForReaction(
  Libnucnet__Net *self,
  Libnucnet__Reaction  *p_reaction,
  double d_t9,
  double d_rho,
  void *p_user_data,
  double *p_forward,
  double *p_reverse
)
{

  /*========================================================================
  // Check input.
  //======================================================================*/

  if( !self ) {
     LIBNUCNET__ERROR( "Invalid network" );
  }

  if( !p_reaction ) {
     LIBNUCNET__ERROR( "Invalid reaction" );
  }

  if( d_t9 < 0. ) {
     LIBNUCNET__ERROR( "Invalid temperature" );
  }

  if( d_rho < 0. ) {
     LIBNUCNET__ERROR( "Invalid density" );
  }

  /*========================================================================
  // Calculate forward rate.
  //======================================================================*/

  *p_forward =
     Libnucnet__Reaction__computeRate(
       p_reaction,
       d_t9,
       p_user_data
     );

  /*========================================================================
  // Get reverse rate.
  //======================================================================*/

  *p_reverse = 
     Libnucnet__Net__computeReverseRate( 
       self, p_reaction, p_forward, d_t9, d_rho
     );

}

/*##############################################################################
// Libnucnet__Net__computeReverseRate().
//############################################################################*/

double
Libnucnet__Net__computeReverseRate(
  Libnucnet__Net *self,
  Libnucnet__Reaction *p_reaction,
  double *p_forward,
  double d_t9,
  double d_rho
) {

  typedef struct {
    Libnucnet__Net *pNet;
    double dT9;
    double dRho;
    double dExp;
  } extra_data;

  double d_exp;
  size_t i_number_nuclide_reactants, i_number_nuclide_products;
  size_t i_number_other_reactants;
  extra_data *p_extra_data;

  /*============================================================================
  // Get number of reactants and products.
  //==========================================================================*/

  i_number_nuclide_reactants =
    (size_t) xmlListSize( p_reaction->pReactantList );

  i_number_nuclide_products =
    (size_t) xmlListSize( p_reaction->pProductList );

  i_number_other_reactants =
    (size_t) xmlListSize( p_reaction->pOtherReactantList );

  /*============================================================================
  // Return 0 if weak reaction or other decay.
  //==========================================================================*/

  if(
    Libnucnet__Reaction__isWeak( p_reaction ) ||
    ( i_number_nuclide_reactants + i_number_other_reactants == 1 )
  )
    return 0;

  /*============================================================================
  // Allocate extra data structure.
  //==========================================================================*/

  p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_extra_data ) LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_extra_data->pNet = self;
  p_extra_data->dT9 = d_t9;
  p_extra_data->dRho = d_rho;
  p_extra_data->dExp = 0.;

  /*============================================================================
  // Loop over reactants and products.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pReactantList,
    (xmlListWalker) Libnucnet__Net__computeReverseRateWalker,
    p_extra_data
  );

  p_extra_data->dExp *= -1.;

  xmlListWalk(
    p_reaction->pProductList,
    (xmlListWalker) Libnucnet__Net__computeReverseRateWalker,
    p_extra_data
  );

  /*============================================================================
  // Get density term.
  //==========================================================================*/

  p_extra_data->dExp +=
    ( 
       ( double ) i_number_nuclide_reactants -
       ( double ) i_number_nuclide_products
    )
    * log( d_rho );

  /*============================================================================
  // Assign to d_exp.
  //==========================================================================*/

  d_exp = p_extra_data->dExp;

  free( p_extra_data );

  /*============================================================================
  // Check factor.
  //==========================================================================*/

  if( d_exp > D_LARGE ) {

    *p_forward = 0.;
    return 0.;

  }

  /*============================================================================
  // Return result.
  //==========================================================================*/

  if( d_exp < -300. ) {

    return 0.;

  } else {

    return *p_forward * exp( d_exp ) *
      Libnucnet__Reaction__getDuplicateProductFactor( p_reaction ) / 
      Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

  }

}

/*##############################################################################
// Libnucnet__Net__computeReverseRateWalker().
//############################################################################*/

int
Libnucnet__Net__computeReverseRateWalker(
  Libnucnet__Reaction__Element *p_element, void *p_data
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    double dT9;
    double dRho;
    double dExp;
  } extra_data;

  extra_data *p_extra_data = ( extra_data * ) p_data;

  p_extra_data->dExp -=
    Libnucnet__Nuc__computeSpeciesNseFactor(
      p_extra_data->pNet->pNuc,
      Libnucnet__Nuc__getSpeciesByName(
        p_extra_data->pNet->pNuc,
        (char *) p_element->sxName
      ),
      p_extra_data->dT9,
      p_extra_data->dRho
    );

  return 1;

}

/*##############################################################################
// Libnucnet__Net__computeReactionQValue().
//############################################################################*/

double
Libnucnet__Net__computeReactionQValue(
  Libnucnet__Net *self, Libnucnet__Reaction *p_reaction
) {

  typedef struct {
    Libnucnet__Net *pNet;
    double dQValue;
    int iType;
  } extra_data;

  double d_q_value;
  extra_data *p_extra_data;

  if( !Libnucnet__Net__isValidReaction( self, p_reaction ) ) {
    LIBNUCNET__ERROR( "Invalid reaction" );
  }

  /*============================================================================
  // Allocate memory and assign extra data.
  //==========================================================================*/

  p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_extra_data ) LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_extra_data->pNet = self;
  p_extra_data->dQValue = 0.;

  /*============================================================================
  // Loop on reactants.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pReactantList,
    (xmlListWalker) Libnucnet__Net__computeReactionQValueWalker,
    p_extra_data
  );

  p_extra_data->dQValue *= -1.;

  /*============================================================================
  // Loop on products.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pProductList,
    (xmlListWalker) Libnucnet__Net__computeReactionQValueWalker,
    p_extra_data
  );

  /*============================================================================
  // Check for beta+ and modify Q.
  //==========================================================================*/

  if( Libnucnet__Reaction__isBetaPlus( p_reaction ) )
    p_extra_data->dQValue -=
      2. * WN_MASS_ELECTRON_MEV;

  /*============================================================================
  // Check for positron capture and modify Q.
  //==========================================================================*/

  if( Libnucnet__Reaction__isPositronCapture( p_reaction ) )
    p_extra_data->dQValue +=
      2. * WN_MASS_ELECTRON_MEV;

  /*============================================================================
  // Return Q Value.
  //==========================================================================*/

  d_q_value = p_extra_data->dQValue;

  free( p_extra_data );

  return d_q_value;

}

/*##############################################################################
// Libnucnet__Net__computeReactionQValueWalker().
//############################################################################*/

int
Libnucnet__Net__computeReactionQValueWalker(
  Libnucnet__Reaction__Element *p_element, void *p_data
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    double dQValue;
    int iType;
  } extra_data;

  extra_data *p_extra_data = ( extra_data * ) p_data;

  p_extra_data->dQValue -=
    Libnucnet__Species__getMassExcess(
      Libnucnet__Nuc__getSpeciesByName(
        p_extra_data->pNet->pNuc,
        (char *) p_element->sxName
      )
    );

  return 1;

}

/*##############################################################################
// Libnucnet__Zone__updateRatesForReaction().
//############################################################################*/

void
Libnucnet__Zone__updateRatesForReaction(
  Libnucnet__Zone *self,
  Libnucnet__Reaction *p_reaction,
  double d_forward,
  double d_reverse
) {

  Libnucnet__Rates *p_rates;

  p_rates = (Libnucnet__Rates *) malloc( sizeof( Libnucnet__Rates ) );

  /*============================================================================
  // Assign rates.
  //==========================================================================*/

  p_rates->dForward = d_forward;
  p_rates->dReverse = d_reverse;

  /*============================================================================
  // Assign to hash.
  //==========================================================================*/

  if ( xmlHashUpdateEntry(
         self->pRatesHash,
         p_reaction->sxReaction,
         p_rates,
         (xmlHashDeallocator) Libnucnet__Rates__free
       )
  ) {
     
    LIBNUCNET__ERROR( "Couldn't update entry" );

  }

}

/*##############################################################################
// Libnucnet__Zone__computeJacobian().
//############################################################################*/

WnMatrix *
Libnucnet__Zone__computeJacobian( Libnucnet__Zone *self )
{

  WnMatrix * p_matrix;

  p_matrix = Libnucnet__Zone__computeJacobianMatrix( self );

  WnMatrix__scaleMatrix( p_matrix, -1. );

  return p_matrix;

}

/*##############################################################################
// Libnucnet__Zone__computeJacobianMatrix().
//############################################################################*/

WnMatrix *
Libnucnet__Zone__computeJacobianMatrix( Libnucnet__Zone *self )
{

  typedef struct {
    Libnucnet__Zone *pZone;
    WnMatrix *pMatrix;
    gsl_vector *pAbunds;
  } extra_data;

  WnMatrix *p_matrix;
  extra_data *p_extra_data;
  Libnucnet__NetView * p_view;

  /*============================================================================
  // Get valid view.
  //==========================================================================*/

  p_view = Libnucnet__Zone__getEvolutionNetView( self );

  /*============================================================================
  // Set extra data structure.
  //==========================================================================*/

  p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  p_extra_data->pZone = self;

  p_extra_data->pMatrix =
    WnMatrix__new(
      Libnucnet__Nuc__getNumberOfSpecies( self->pNet->pNuc ), 
      Libnucnet__Nuc__getNumberOfSpecies( self->pNet->pNuc )
    );

  p_extra_data->pAbunds = Libnucnet__Zone__getAbundances( self );

  xmlHashScan(
    p_view->pNet->pReac->pReactionHash,
    (xmlHashScanner) Libnucnet__Zone__computeJacobianMatrixCallback,
    p_extra_data
  );

  p_matrix = p_extra_data->pMatrix;

  gsl_vector_free( p_extra_data->pAbunds );
  free( p_extra_data );

  return p_matrix;

}
 
/*##############################################################################
// Libnucnet__Zone__computeJacobianMatrixCallback().
//############################################################################*/

void
Libnucnet__Zone__computeJacobianMatrixCallback(
  Libnucnet__Reaction *p_reaction,
  void *p_data,
  xmlChar *sx_reaction
) {

  typedef struct {
    Libnucnet__Zone *pZone;
    WnMatrix *pMatrix;
    gsl_vector *pAbunds;
  } extra_data;

  typedef struct 
  {
    Libnucnet__Species **pSpecies;
    Libnucnet__Nuc * pNuc;
    int iCount;
  } species_data;

  species_data * p_reactants, * p_products;
  Libnucnet__Rates *p_rates;

  extra_data *p_extra_data = ( extra_data * ) p_data;
  size_t i, j, i_number_reactants, i_number_products;
  double d_value;

  /*============================================================================
  // Get rates.
  //==========================================================================*/

  if( 
      !(
         p_rates =
           (Libnucnet__Rates * )
            xmlHashLookup( p_extra_data->pZone->pRatesHash, sx_reaction )
      )
  ) {
     LIBNUCNET__ERROR( "No rates in structure for Jacobian" );
  }

  /*============================================================================
  // Create and populate arrays.
  //==========================================================================*/

  p_reactants = ( species_data * ) malloc( sizeof( species_data ) );
  p_products = ( species_data * ) malloc( sizeof( species_data ) );

  if( !p_reactants || !p_products )
    LIBNUCNET__ERROR( "Couldn't allocate memory for species array" );

  p_reactants->pNuc = p_extra_data->pZone->pNet->pNuc;
  p_products->pNuc = p_extra_data->pZone->pNet->pNuc;

  i_number_reactants = (size_t) xmlListSize( p_reaction->pReactantList );
  i_number_products = (size_t) xmlListSize( p_reaction->pProductList );

  p_reactants->pSpecies =
    (Libnucnet__Species **)
    malloc( sizeof( Libnucnet__Species * ) * i_number_reactants );

  p_products->pSpecies =
    (Libnucnet__Species **)
    malloc( sizeof( Libnucnet__Species * ) * i_number_products );

  if( !p_reactants->pSpecies || !p_products->pSpecies )
    LIBNUCNET__ERROR( "Couldn't allocate memory for species array" );

  p_reactants->iCount = 0;
  xmlListWalk(
    p_reaction->pReactantList,
    (xmlListWalker) Libnucnet__Net__populate_reaction_species_array,
    p_reactants
  );

  p_products->iCount = 0;
  xmlListWalk(
    p_reaction->pProductList,
    (xmlListWalker) Libnucnet__Net__populate_reaction_species_array,
    p_products
  );

  /*============================================================================
  // Loop on reactants.
  //==========================================================================*/

  for( i = 0; i < i_number_reactants; i++ )
  {

    d_value =
      p_rates->dForward /
        Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

    for( j = 0; j < i_number_reactants; j++ )
    {

      if( i != j )
      {

        d_value *=
          gsl_vector_get(
            p_extra_data->pAbunds,
            p_reactants->pSpecies[j]->iIndex
          );

      }

    }

    for( j = 0; j < i_number_reactants; j++ )
    {

      WnMatrix__assignElement(
	p_extra_data->pMatrix,
	Libnucnet__Species__getIndex( p_reactants->pSpecies[j] ) + 1,
	Libnucnet__Species__getIndex( p_reactants->pSpecies[i] ) + 1,
	d_value
      );

    }

    for( j = 0; j < i_number_products; j++ )
    {

      WnMatrix__assignElement(
	p_extra_data->pMatrix,
	Libnucnet__Species__getIndex( p_products->pSpecies[j] ) + 1,
	Libnucnet__Species__getIndex( p_reactants->pSpecies[i] ) + 1,
	-d_value
      );

    }

  }

  /*============================================================================
  // Loop on products.
  //==========================================================================*/

  for( i = 0; i < i_number_products; i++ )
  {

    d_value =
      p_rates->dReverse /
        Libnucnet__Reaction__getDuplicateProductFactor( p_reaction );

    for( j = 0; j < i_number_products; j++ )
    {

      if( i != j )
      {

        d_value *=
          gsl_vector_get(
            p_extra_data->pAbunds,
            p_products->pSpecies[j]->iIndex
          );

      }

    }

    for( j = 0; j < i_number_reactants; j++ )
    {

      WnMatrix__assignElement(
	p_extra_data->pMatrix,
	Libnucnet__Species__getIndex( p_reactants->pSpecies[j] ) + 1,
	Libnucnet__Species__getIndex( p_products->pSpecies[i] ) + 1,
	-d_value
      );

    }

    for( j = 0; j < i_number_products; j++ )
    {

      WnMatrix__assignElement(
	p_extra_data->pMatrix,
	Libnucnet__Species__getIndex( p_products->pSpecies[j] ) + 1,
	Libnucnet__Species__getIndex( p_products->pSpecies[i] ) + 1,
	d_value
      );

    }

  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  free( p_reactants->pSpecies );
  free( p_reactants );

  free( p_products->pSpecies );
  free( p_products );

}

/*##############################################################################
// Libnucnet__Zone__computeFlowVector
//############################################################################*/

gsl_vector *
Libnucnet__Zone__computeFlowVector( Libnucnet__Zone *self ) {

  typedef struct {
    Libnucnet__Zone *pZone;
    gsl_vector *pAbunds;
    gsl_vector *pFlow;
  } extra_data;

  gsl_vector *p_flow;
  Libnucnet__NetView * p_view;

  extra_data * p_extra_data;

  /*============================================================================
  // Get valid view.
  //==========================================================================*/

  p_view = Libnucnet__Zone__getEvolutionNetView( self );

  /*============================================================================
  // Set extra data structure.
  //==========================================================================*/

  p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  p_extra_data->pZone = self;

  p_extra_data->pAbunds = Libnucnet__Zone__getAbundances( self );

  p_extra_data->pFlow =
    gsl_vector_calloc( 
      Libnucnet__Nuc__getNumberOfSpecies( self->pNet->pNuc )
    );

  /*============================================================================
  // Hash callback.
  //==========================================================================*/

  xmlHashScan(
    p_view->pNet->pReac->pReactionHash,
    (xmlHashScanner) Libnucnet__Zone__computeFlowVectorCallback,
    p_extra_data
  );

  p_flow = p_extra_data->pFlow;

  gsl_vector_free( p_extra_data->pAbunds );
  free( p_extra_data );

  return p_flow;

}

/*##############################################################################
// Libnucnet__Zone__computeFlowVectorCallback()
//############################################################################*/

void
Libnucnet__Zone__computeFlowVectorCallback(
  Libnucnet__Reaction *p_reaction,
  void *p_data,
  xmlChar *sx_reaction
) {

  typedef struct {
    Libnucnet__Zone *pZone;
    gsl_vector *pAbunds;
    gsl_vector *pFlow;
  } extra_data;

  typedef struct 
  {
    Libnucnet__Species **pSpecies;
    Libnucnet__Nuc *pNuc;
    int iCount;
  } species_data;

  double d_u;
  Libnucnet__Rates * p_rates;
  species_data * p_reactants, * p_products;
  size_t i, i_number_reactants, i_number_products;
  double d_value;

  extra_data * p_extra_data = ( extra_data * ) p_data;

  /*============================================================================
  // Get rates for reaction.
  //==========================================================================*/

  if( 
      !(
         p_rates =
           (Libnucnet__Rates * )
            xmlHashLookup( p_extra_data->pZone->pRatesHash, sx_reaction )
      )
  ) {
     LIBNUCNET__ERROR( "No rates in structure for flow vector" );
  }

  /*============================================================================
  // Create and populate arrays.
  //==========================================================================*/

  p_reactants = ( species_data * ) malloc( sizeof( species_data ) );
  p_products = ( species_data * ) malloc( sizeof( species_data ) );

  if( !p_reactants || !p_products )
    LIBNUCNET__ERROR( "Couldn't allocate memory for species array" );

  p_reactants->pNuc = p_extra_data->pZone->pNet->pNuc;
  p_products->pNuc = p_extra_data->pZone->pNet->pNuc;

  i_number_reactants = (size_t) xmlListSize( p_reaction->pReactantList );
  i_number_products = (size_t) xmlListSize( p_reaction->pProductList );

  p_reactants->pSpecies =
    (Libnucnet__Species **)
    malloc( sizeof( Libnucnet__Species * ) * i_number_reactants );

  p_products->pSpecies =
    (Libnucnet__Species **)
    malloc( sizeof( Libnucnet__Species * ) * i_number_products );

  if( !p_reactants->pSpecies || !p_products->pSpecies )
    LIBNUCNET__ERROR( "Couldn't allocate memory for species array" );

  p_reactants->iCount = 0;
  xmlListWalk(
    p_reaction->pReactantList,
    (xmlListWalker) Libnucnet__Net__populate_reaction_species_array,
    p_reactants
  );

  p_products->iCount = 0;
  xmlListWalk(
    p_reaction->pProductList,
    (xmlListWalker) Libnucnet__Net__populate_reaction_species_array,
    p_products
  );

  /*============================================================================
  // First forward part of reaction.
  //==========================================================================*/

  d_value =
    p_rates->dForward /
      Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

  for( i = 0; i < i_number_reactants; i++ )
  {

    d_value *=
      gsl_vector_get(
        p_extra_data->pAbunds,
        p_reactants->pSpecies[i]->iIndex
      );

  }

  d_u = d_value;

  /*============================================================================
  // Now reverse part of reaction.
  //==========================================================================*/

  d_value =
    p_rates->dReverse /
      Libnucnet__Reaction__getDuplicateProductFactor( p_reaction );

  for( i = 0; i < i_number_products; i++ )
  {

    d_value *=
      gsl_vector_get(
        p_extra_data->pAbunds,
        p_products->pSpecies[i]->iIndex
      );

  }

  d_u -= d_value;

  /*============================================================================
  // Now assign to vector.
  //==========================================================================*/

  for( i = 0; i < i_number_reactants; i++ )
  {
    p_extra_data->pFlow->data[p_reactants->pSpecies[i]->iIndex] -= d_u;
  }

  for( i = 0; i < i_number_products; i++ )
  {
    p_extra_data->pFlow->data[p_products->pSpecies[i]->iIndex] += d_u;
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  free( p_reactants->pSpecies );
  free( p_reactants );

  free( p_products->pSpecies );
  free( p_products );

}

/*##############################################################################
// Libnucnet__Zone__updateTimeStep().
//############################################################################*/

void
Libnucnet__Zone__updateTimeStep(
  Libnucnet__Zone *self,
  double *p_dt,
  double d_regt,
  double d_regy,
  double d_ymin
)
{

  typedef struct {
    const Libnucnet__Zone *pZone;
    double dRegy;
    double dYmin;
    double dDt;
    double dDtnew;
  } extra_data;

  extra_data *p_extra_data = ( extra_data * ) malloc( sizeof( extra_data) );

  p_extra_data->dDt = *p_dt;

  p_extra_data->pZone = self;
  p_extra_data->dRegy = d_regy;
  p_extra_data->dYmin = d_ymin;
  p_extra_data->dDtnew = ( 1. + d_regt ) * p_extra_data->dDt;

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__updateTimeStepCallback,
    p_extra_data
  );

  *p_dt = p_extra_data->dDtnew;

  free( p_extra_data );

}

/*##############################################################################
// Libnucnet__Zone__updateTimeStepCallback().
//############################################################################*/

void
Libnucnet__Zone__updateTimeStepCallback(
  Libnucnet__Species *p_species,
  void *p_data,
  xmlChar *sx_species
)
{

  typedef struct {
    const Libnucnet__Zone *pZone;
    double dRegy;
    double dYmin;
    double dDt;
    double dDtnew;
  } extra_data;

  double d_y, d_dy, d_check;
  double d_tiny = 1.e-300;

  extra_data *p_extra_data = ( extra_data * ) p_data;

  if( !sx_species )
    LIBNUCNET__ERROR( "No such species" );

  d_y = Libnucnet__Zone__getSpeciesAbundance( p_extra_data->pZone, p_species );
  d_dy =
    Libnucnet__Zone__getSpeciesAbundanceChange(
      p_extra_data->pZone, p_species
    );

  if ( d_y > p_extra_data->dYmin && !WnMatrix__value_is_zero( d_dy ) ) {

    d_check =
      p_extra_data->dDt *
      p_extra_data->dRegy *
      fabs( d_y / ( d_dy + d_tiny) );

    if ( d_check < p_extra_data->dDtnew ) {
      p_extra_data->dDtnew = d_check;
    }

  }

}

/*##############################################################################
// Libnucnet__getNumberOfZones().
//############################################################################*/

size_t Libnucnet__getNumberOfZones( const Libnucnet *self ) {

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid Network" );
  }

  return (size_t) xmlHashSize( self->pZoneHash );

}

/*##############################################################################
// Libnucnet__removeZone()
//############################################################################*/

void
Libnucnet__removeZone(
  Libnucnet *self, Libnucnet__Zone * p_zone
){

  xmlHashRemoveEntry3(
    self->pZoneHash,
    p_zone->sxLabel[0],
    p_zone->sxLabel[1],
    p_zone->sxLabel[2],
    ( xmlHashDeallocator )Libnucnet__Zone__hashFree
  );

}

/*##############################################################################
// Libnucnet__Zone__updateSpeciesAbundance().
//############################################################################*/

void
Libnucnet__Zone__updateSpeciesAbundance(
  Libnucnet__Zone *self,
  Libnucnet__Species *p_species,
  double d_y
) {

  double *p_abund_new;

  if( !self || !p_species ) LIBNUCNET__ERROR( "Invalid input" );

  if( !WnMatrix__value_is_zero( d_y ) ) {

    p_abund_new =
      ( double * ) malloc( sizeof( double ) );

    *p_abund_new = d_y;

    if(
        xmlHashUpdateEntry(
           self->pAbundanceHash,
           p_species->sxName,
           p_abund_new,
           (xmlHashDeallocator) Libnucnet__Abundance__free
        ) != 0
    ) {
        LIBNUCNET__ERROR( "Couldn't update abundance" );
    }

  } else {
    
    if(
        xmlHashLookup(
          self->pAbundanceHash,
          p_species->sxName
        )
    ) {
        xmlHashRemoveEntry(
          self->pAbundanceHash,
          p_species->sxName,
          (xmlHashDeallocator) Libnucnet__Abundance__free
        );
    }

  }

}

/*##############################################################################
// Libnucnet__Zone__updateSpeciesAbundanceChange().
//############################################################################*/

void
Libnucnet__Zone__updateSpeciesAbundanceChange(
  Libnucnet__Zone *self, Libnucnet__Species *p_species, double d_dy
) {

  double *p_abund_change_new;

  if( !self || !p_species ) LIBNUCNET__ERROR( "Invalid input" );

  if( !WnMatrix__value_is_zero( d_dy ) ) {

    p_abund_change_new =
      ( double * ) malloc( sizeof( double ) );

    *p_abund_change_new = d_dy;

    if(
        xmlHashUpdateEntry(
           self->pAbundanceChangeHash,
           p_species->sxName,
           p_abund_change_new,
           (xmlHashDeallocator) Libnucnet__AbundanceChange__free
        ) != 0
    ) {
        LIBNUCNET__ERROR( "Couldn't update abundance change" );
    }

  } else {

    if(
        xmlHashLookup(
           self->pAbundanceChangeHash,
           p_species->sxName
        )
    ) {
        xmlHashRemoveEntry(
           self->pAbundanceChangeHash,
           p_species->sxName,
           (xmlHashDeallocator) Libnucnet__AbundanceChange__free
        );
    }

  }
        

}

/*##############################################################################
// Libnucnet__Zone__getSpeciesAbundance().
//############################################################################*/

double
Libnucnet__Zone__getSpeciesAbundance(
  const Libnucnet__Zone *self,
  const Libnucnet__Species *p_species
) {

  double *p_abund;

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  if( !p_species ) LIBNUCNET__ERROR( "Invalid species" );

  p_abund =
     ( double * )
     xmlHashLookup(
         self->pAbundanceHash,
         p_species->sxName
     );

  if ( !p_abund ) {

    return 0.;

  }

  return *p_abund;

}

/*##############################################################################
// Libnucnet__Zone__getSpeciesAbundanceChange().
//############################################################################*/

double
Libnucnet__Zone__getSpeciesAbundanceChange(
  const Libnucnet__Zone *self,
  const Libnucnet__Species *p_species
) {

  double *p_abund_change;

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  if( !p_species ) LIBNUCNET__ERROR( "Invalid species" );

  p_abund_change =
     ( double * )
     xmlHashLookup(
         self->pAbundanceChangeHash,
         p_species->sxName
     );

  if ( !p_abund_change ) {

    return 0.;

  }

  return *p_abund_change;

}

/*##############################################################################
// Libnucnet__getZoneByLabels().
//############################################################################*/

Libnucnet__Zone *
Libnucnet__getZoneByLabels(
  const Libnucnet *self,
  const char *s_label_1,
  const char *s_label_2,
  const char *s_label_3
) {

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid structure" );
  }

  return
    ( Libnucnet__Zone * )
    xmlHashLookup3(
      self->pZoneHash, 
      (const xmlChar *) s_label_1, 
      (const xmlChar *) s_label_2, 
      (const xmlChar *)  s_label_3 
    );

}

/*##############################################################################
// Libnucnet__Zone__free()
//############################################################################*/

void
Libnucnet__Zone__free( Libnucnet__Zone *self ) {

  int i;

  xmlHashFree(
    self->pAbundanceHash,
    (xmlHashDeallocator) Libnucnet__Abundance__free
  );
  xmlHashFree(
    self->pAbundanceChangeHash,
    (xmlHashDeallocator) Libnucnet__AbundanceChange__free
  );
  xmlHashFree( self->pRatesHash, (xmlHashDeallocator) Libnucnet__Rates__free );

  xmlHashFree( self->pPropertyHash, (xmlHashDeallocator) xmlFree );

  xmlHashFree(
    self->pViewHash,
    (xmlHashDeallocator) Libnucnet__NetView__free
  );

  xmlHashScan(
    self->pNet->pReac->pFDHash,
    (xmlHashScanner) Libnucnet__Zone__user_rate_function_data_free,
    self
  );

  xmlHashFree( self->pUserDataHash, NULL );

  for( i = 0; i < MAX_ZONE_LABELS; i++ ) {
    free( self->sxLabel[i] );
  }

  free( self );
  
  self = NULL;

  return;

}

/*##############################################################################
// Libnucnet__Zone__user_data_free()
//############################################################################*/

void
Libnucnet__Zone__user_rate_function_data_free(
  Libnucnet__Reac__FD *p_fd,
  Libnucnet__Zone *p_zone,
  xmlChar *sx_function_key
)
{

  void *p_data;

  if( p_fd->pfDeallocator )
  {

    p_data =
      Libnucnet__Zone__getDataForUserRateFunction(
        p_zone,
        (const char *) sx_function_key
      );

    if( p_data ) p_fd->pfDeallocator( p_data );

  } 

}

/*##############################################################################
// Libnucnet__Zone__hashFree()
//############################################################################*/

void
Libnucnet__Zone__hashFree( Libnucnet__Zone *self, xmlChar *sx_label ) {

  if( !sx_label ) {
    LIBNUCNET__ERROR( "Trying to free invalid zone" );
  }

  Libnucnet__Zone__free( self );

}

/*##############################################################################
// Libnucnet__Abundance__free()
//############################################################################*/

void Libnucnet__Abundance__free( double *self, xmlChar *sx_name ) {

  if( !sx_name ) {
      LIBNUCNET__ERROR( "Trying to free non-existent species" );
  }

  free( self );

  return;

}

/*##############################################################################
// Libnucnet__AbundanceChange__free()
//############################################################################*/

void Libnucnet__AbundanceChange__free( double *self, xmlChar *s_name ) {

  if( !s_name ) {
      LIBNUCNET__ERROR( "Trying to free non-existent species" );
  }

  free( self );

  return;

}

/*##############################################################################
// Libnucnet__Net__isValidReaction()
//############################################################################*/

int
Libnucnet__Net__isValidReaction(
  const Libnucnet__Net *self,
  const Libnucnet__Reaction *p_reaction
)
{

  typedef struct {
    const Libnucnet__Net *pNet;
    int iValid;
  } extra_data;

  extra_data *p_extra_data;

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid network" );
  }

  if( !p_reaction ) {
    LIBNUCNET__ERROR( "Invalid reaction" );
  }

  /*============================================================================
  // Allocate extra data structure.
  //==========================================================================*/

  p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_extra_data ) LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_extra_data->pNet = self;
  p_extra_data->iValid = 1;

  /*============================================================================
  // Check reactants.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pReactantList,
    (xmlListWalker) Libnucnet__Net__isValidReactionWalker,
    p_extra_data
  );

  if( !p_extra_data->iValid ) {
    free( p_extra_data );
    return 0;
  }

  /*============================================================================
  // Check products.
  //==========================================================================*/

  xmlListWalk(
    p_reaction->pProductList,
    (xmlListWalker) Libnucnet__Net__isValidReactionWalker,
    p_extra_data
  );

  if( !p_extra_data->iValid ) {
    free( p_extra_data );
    return 0;
  }

  free( p_extra_data );

  /*----------------------------------------------------------------------
  // For valid reaction, check charge, baryon, and lepton conservation.
  //--------------------------------------------------------------------*/

  Libnucnet__Net__checkReactionForBaryonConservation(
    self,
    p_reaction
  );
  Libnucnet__Net__checkReactionForChargeConservation(
    self,
    p_reaction
  );
  Libnucnet__Net__checkReactionForLeptonConservation(
    self,
    p_reaction
  );

  return 1;

}

/*##############################################################################
// Libnucnet__Net__isValidReactionWalker()
//############################################################################*/

int
Libnucnet__Net__isValidReactionWalker(
  Libnucnet__Reaction__Element *p_element, void *p_data
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    int iValid;
  } extra_data;

  extra_data *p_extra_data = ( extra_data * ) p_data;

  if(
    !( Libnucnet__Nuc__getSpeciesByName(
         p_extra_data->pNet->pNuc,
         (char *) p_element->sxName
       )
    )
  ) { 
    p_extra_data->iValid = 0;
    return 0;
  }

  return 1;

}

/*##############################################################################
// Libnucnet__freeAllZones()
//############################################################################*/

void
Libnucnet__freeAllZones( Libnucnet *self )
{

  if( !self ) LIBNUCNET__ERROR( "Invalid input" );

  xmlHashFree(
    self->pZoneHash,
    (xmlHashDeallocator) Libnucnet__Zone__hashFree
  );

  Libnucnet__initializeHash( self );

}

/*##############################################################################
// Libnucnet__Zone__getLabel()
//############################################################################*/

const char *
Libnucnet__Zone__getLabel( const Libnucnet__Zone *self, int i_label ) {

  if( !self ) LIBNUCNET__ERROR( "Trying to get label of invalid zone" );

  if( 
      (
         i_label < 1 || i_label > MAX_ZONE_LABELS ) 
      )
  {
    LIBNUCNET__ERROR( "Invalid label number" );
  }

  return (const char *) self->sxLabel[i_label-1];

}

/*##############################################################################
// Libnucnet__Rates__free()
//############################################################################*/

void Libnucnet__Rates__free( Libnucnet__Rates *self, xmlChar *sx_reaction ) {

  if( !sx_reaction ) {
    LIBNUCNET__ERROR( "Trying to free invalid reaction" );
  }

  free( self );

}

/*##############################################################################
// Libnucnet__Net__free()
//############################################################################*/

void
Libnucnet__Net__free( Libnucnet__Net *self ) {

  if( !self ) return;

  Libnucnet__Nuc__free( self->pNuc );
  Libnucnet__Reac__free( self->pReac );

  free( self );

}

/*##############################################################################
// Libnucnet__NetView__free()
//############################################################################*/

void
Libnucnet__NetView__free( Libnucnet__NetView *self ) {

  if( !self ) return;

  Libnucnet__Net__free( self->pNet );

  free( self );

}

/*##############################################################################
// Libnucnet__Zone__getNetView()
//############################################################################*/

Libnucnet__NetView *
Libnucnet__Zone__getNetView(
  Libnucnet__Zone * self,
  const char * s_label1,
  const char * s_label2,
  const char * s_label3
)
{

  Libnucnet__NetView * p_view;

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  /*============================================================================
  // Lookup view and return.
  //==========================================================================*/

  p_view =
    (Libnucnet__NetView *)
    xmlHashLookup3(
      self->pViewHash,
      (const xmlChar *) s_label1,
      (const xmlChar *) s_label2,
      (const xmlChar *) s_label3
    );

  return p_view;

}

/*##############################################################################
// Libnucnet__NetView__wasNetUpdated()
//############################################################################*/

int
Libnucnet__NetView__wasNetUpdated(
  Libnucnet__NetView * self
)
{

  if( !self ) LIBNUCNET__ERROR( "Invalid net view" );

  if(
      self->pParentNet->pNuc->iUpdate == self->pNet->pNuc->iUpdate &&
      self->pParentNet->pReac->iUpdate == self->pNet->pReac->iUpdate
  )
    return 0;
  else
    return 1;

}

/*##############################################################################
// Libnucnet__Zone__updateNetView().
//############################################################################*/

int
Libnucnet__Zone__updateNetView(
  Libnucnet__Zone * self,
  const char * s_label1,
  const char * s_label2,
  const char * s_label3,
  Libnucnet__NetView * p_view
)
{

  if( !self ) LIBNUCNET__ERROR( "Updating net view of invalid zone" );

  return
    xmlHashUpdateEntry3(
      self->pViewHash,
      (const xmlChar *) s_label1,
      (const xmlChar *) s_label2,
      (const xmlChar *) s_label3,
      p_view,
      (xmlHashDeallocator) Libnucnet__NetView__free
    ) + 1;

}

/*##############################################################################
// Libnucnet__Zone__removeNetView().
//############################################################################*/

int
Libnucnet__Zone__removeNetView(
  Libnucnet__Zone * self,
  const char * s_label1,
  const char * s_label2,
  const char * s_label3
)
{

  return
    xmlHashRemoveEntry3(
      self->pViewHash,
      (const xmlChar *) s_label1,
      (const xmlChar *) s_label2,
      (const xmlChar *) s_label3,
      (xmlHashDeallocator) Libnucnet__NetView__free
    ) + 1;

}

/*##############################################################################
// Libnucnet__NetView__addReaction()
//############################################################################*/

int
Libnucnet__NetView__addReaction(
  Libnucnet__NetView * p_net_view,
  Libnucnet__Reaction * p_reaction
)
{

  if( !p_net_view || !p_reaction ) LIBNUCNET__ERROR( "Invalid input" );

  if(
    Libnucnet__Net__isValidReaction(
      Libnucnet__NetView__getNet( p_net_view ),
      p_reaction
    )
  )
  {
    return
      xmlHashAddEntry(
        p_net_view->pNet->pReac->pReactionHash,
        (const xmlChar *) Libnucnet__Reaction__getString( p_reaction ),
        p_reaction
      ) + 1;
  }
  else
  {
    return 0;
  }

}

/*##############################################################################
// Libnucnet__NetView__removeReaction()
//############################################################################*/

int
Libnucnet__NetView__removeReaction(
  Libnucnet__NetView * p_net_view,
  Libnucnet__Reaction * p_reaction
)
{

  if( !p_net_view || !p_reaction ) LIBNUCNET__ERROR( "Invalid input" );

  return
    xmlHashRemoveEntry(
      p_net_view->pNet->pReac->pReactionHash,
      (const xmlChar *) Libnucnet__Reaction__getString( p_reaction ),
      NULL
    ) + 1;

}

/*##############################################################################
// Libnucnet__free()
//############################################################################*/

void
Libnucnet__free( Libnucnet *self ) {

  if( !self ) return;

  xmlHashFree(
    self->pZoneHash,
    (xmlHashDeallocator) Libnucnet__Zone__hashFree
  );

  Libnucnet__Net__free( self->pNet );

  xmlFree( self->sxMassFractionFormat );
  
  free( self );

}

/*##############################################################################
// Libnucnet__Net__getReac()
//############################################################################*/

Libnucnet__Reac * Libnucnet__Net__getReac( const Libnucnet__Net *self )
{

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid network" );
  }

  return self->pReac;

}

/*##############################################################################
// Libnucnet__Net__getNuc()
//############################################################################*/

Libnucnet__Nuc * Libnucnet__Net__getNuc( const Libnucnet__Net *self )
{

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid network" );
  }

  return self->pNuc;

}

/*##############################################################################
// Libnucnet__Zone__getNet()
//############################################################################*/

Libnucnet__Net * Libnucnet__Zone__getNet( const Libnucnet__Zone *self )
{

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid zone" );
  }

  if( !self->pNet ) {
    LIBNUCNET__ERROR( "Invalid network" );
  }

  return self->pNet;

}

/*##############################################################################
// Libnucnet__getNet()
//############################################################################*/

Libnucnet__Net *Libnucnet__getNet( Libnucnet *self )
{

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid Libnucnet structure" );
  }

  return self->pNet;

}

/*##############################################################################
// Libnucnet__Net__is_valid_input_xml()
//############################################################################*/

int Libnucnet__Net__is_valid_input_xml( const char *s_filename ) {

  xmlDocPtr p_doc;
  xmlSchemaParserCtxtPtr p_parser_ctxt;
  xmlSchemaValidCtxtPtr p_valid_ctxt;
  xmlSchemaPtr p_schema;
  int i_valid;

  p_parser_ctxt = xmlSchemaNewParserCtxt( LIBNUCNET__NET__SCHEMA );

  p_schema = xmlSchemaParse( p_parser_ctxt );

  if( !p_schema ) {
     LIBNUCNET__ERROR( "Schema not found or invalid" );
  }

  p_valid_ctxt = xmlSchemaNewValidCtxt( p_schema );

  p_doc = xmlParseFile( s_filename );

  if( xmlXIncludeProcess( p_doc ) == -1 )
    LIBNUCNET__ERROR( "Problem including xml" );

  if( xmlSchemaValidateDoc( p_valid_ctxt, p_doc ) == 0 ) {
     i_valid = 1;
  } else {
     i_valid = 0;
  }

  xmlFreeDoc( p_doc );
  xmlSchemaFreeValidCtxt( p_valid_ctxt );
  xmlSchemaFree( p_schema );
  xmlSchemaFreeParserCtxt( p_parser_ctxt );

  xmlCleanupParser();

  return i_valid;

}

/*##############################################################################
// Libnucnet__is_valid_input_xml()
//############################################################################*/

int Libnucnet__is_valid_input_xml( const char *s_filename ) {

  xmlDocPtr p_doc;
  xmlSchemaParserCtxtPtr p_parser_ctxt;
  xmlSchemaValidCtxtPtr p_valid_ctxt;
  xmlSchemaPtr p_schema;
  int i_valid;

  p_parser_ctxt = xmlSchemaNewParserCtxt( LIBNUCNET__SCHEMA );

  p_schema = xmlSchemaParse( p_parser_ctxt );

  if( !p_schema ) {
     LIBNUCNET__ERROR( "Schema not found or invalid" );
  }

  p_valid_ctxt = xmlSchemaNewValidCtxt( p_schema );

  p_doc = xmlParseFile( s_filename );

  if( xmlXIncludeProcess( p_doc ) == -1 )
    LIBNUCNET__ERROR( "Problem including xml" );

  if( xmlSchemaValidateDoc( p_valid_ctxt, p_doc ) == 0 ) {
     i_valid = 1;
  } else {
     i_valid = 0;
  }

  xmlFreeDoc( p_doc );
  xmlSchemaFreeValidCtxt( p_valid_ctxt );
  xmlSchemaFree( p_schema );
  xmlSchemaFreeParserCtxt( p_parser_ctxt );

  xmlCleanupParser();

  return i_valid;

}

/*##############################################################################
// Libnucnet__is_valid_zone_data_xml()
//############################################################################*/

int Libnucnet__is_valid_zone_data_xml( const char *s_filename ) {

  xmlDocPtr p_doc;
  xmlSchemaParserCtxtPtr p_parser_ctxt;
  xmlSchemaValidCtxtPtr p_valid_ctxt;
  xmlSchemaPtr p_schema;
  int i_valid;

  p_parser_ctxt = xmlSchemaNewParserCtxt( LIBNUCNET__ZONE_DATA__SCHEMA );

  p_schema = xmlSchemaParse( p_parser_ctxt );

  if( !p_schema ) {
     LIBNUCNET__ERROR( "Schema not found or invalid" );
  }

  p_valid_ctxt = xmlSchemaNewValidCtxt( p_schema );

  p_doc = xmlParseFile( s_filename );

  if( xmlXIncludeProcess( p_doc ) == -1 )
    LIBNUCNET__ERROR( "Problem including xml" );

  if( xmlSchemaValidateDoc( p_valid_ctxt, p_doc ) == 0 ) {
     i_valid = 1;
  } else {
     i_valid = 0;
  }

  xmlFreeDoc( p_doc );
  xmlSchemaFreeValidCtxt( p_valid_ctxt );
  xmlSchemaFree( p_schema );
  xmlSchemaFreeParserCtxt( p_parser_ctxt );

  xmlCleanupParser();

  return i_valid;

}

/*##############################################################################
// Libnucnet__Zone__computeAMoment()
//############################################################################*/

double
Libnucnet__Zone__computeAMoment( const Libnucnet__Zone *self, int n )
{
  
  typedef struct {
    double dResult;
    int n;
    const Libnucnet__Zone *pZone;
  } extra_data;

  double d_result;
  extra_data *p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_extra_data )
    LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_extra_data->dResult = 0.;
  p_extra_data->n = n;
  p_extra_data->pZone = self;

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid zone" );
  }

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__computeAMomentCallback,
    p_extra_data
  );
    
  d_result = p_extra_data->dResult;
  free( p_extra_data );

  return d_result;

}

/*##############################################################################
// Libnucnet__Zone__computeAMomentCallback()
//############################################################################*/

void
Libnucnet__Zone__computeAMomentCallback(
  Libnucnet__Species *p_species, void *p_data, xmlChar *sx_species 
)
{

  typedef struct {
    double dResult;
    int n;
    const Libnucnet__Zone *pZone;
  } extra_data;

  extra_data *p_extra_data = ( extra_data * ) p_data;

  if( !sx_species )
    LIBNUCNET__ERROR( "No such species" );

  p_extra_data->dResult +=
    pow(
         (double) Libnucnet__Species__getA( p_species ),
         (double) p_extra_data->n
    ) * Libnucnet__Zone__getSpeciesAbundance( p_extra_data->pZone, p_species );

} 

/*##############################################################################
// Libnucnet__Zone__computeZMoment()
//############################################################################*/

double
Libnucnet__Zone__computeZMoment( const Libnucnet__Zone *self, int n )
{
  
  typedef struct {
    double dResult;
    int n;
    const Libnucnet__Zone *pZone;
  } extra_data;

  double d_result;
  extra_data *p_extra_data = ( extra_data * ) malloc( sizeof( extra_data ) );

  if( !p_extra_data )
    LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_extra_data->dResult = 0.;
  p_extra_data->n = n;
  p_extra_data->pZone = self;

  if( !self ) {
    LIBNUCNET__ERROR( "Invalid zone" );
  }

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__computeZMomentCallback,
    p_extra_data
  );
    
  d_result = p_extra_data->dResult;
  free( p_extra_data );

  return d_result;

}

/*##############################################################################
// Libnucnet__Zone__computeZMomentCallback()
//############################################################################*/

void
Libnucnet__Zone__computeZMomentCallback(
  Libnucnet__Species *p_species, void *p_data, xmlChar *sx_species 
)
{

  typedef struct {
    double dResult;
    int n;
    const Libnucnet__Zone *pZone;
  } extra_data;

  extra_data *p_extra_data = ( extra_data * ) p_data;

  if( !sx_species )
    LIBNUCNET__ERROR( "No such species" );

  p_extra_data->dResult +=
    pow(
         (double) Libnucnet__Species__getZ( p_species ),
         (double) p_extra_data->n
    ) * Libnucnet__Zone__getSpeciesAbundance( p_extra_data->pZone, p_species );

} 

/*##############################################################################
// Libnucnet__Zone__clearScreeningFunction()
//############################################################################*/

void
Libnucnet__Zone__clearScreeningFunction( Libnucnet__Zone *self )
{
  
  self->pfScreeningFunction = NULL;
  self->pScreeningData = NULL;

}

/*##############################################################################
// Libnucnet__Zone__clearNseCorrectionFactorFunction()
//############################################################################*/

void
Libnucnet__Zone__clearNseCorrectionFactorFunction(
  Libnucnet__Zone *self
)
{
  
  self->pfNseCorrectionFactorFunction = NULL;
  self->pNseCorrectionFactorData = NULL;

}

/*##############################################################################
// Libnucnet__Zone__getRatesForReaction().
//############################################################################*/

void
Libnucnet__Zone__getRatesForReaction(
  const Libnucnet__Zone *self,
  const Libnucnet__Reaction *p_reaction,
  double *p_forward,
  double *p_reverse
) {

  Libnucnet__Rates *p_rates;

  /*============================================================================
  // Retrieve rates from hash.
  //==========================================================================*/

  p_rates =
    ( Libnucnet__Rates * )
    xmlHashLookup( 
      self->pRatesHash, p_reaction->sxReaction 
    );

  /*============================================================================
  // Assign rates.
  //==========================================================================*/

  if( !p_rates )
  {
    fprintf(
      stderr,
      "Reaction: %s\n",
      Libnucnet__Reaction__getString( p_reaction )
    );
    LIBNUCNET__ERROR( "Rates not assigned" );
  }
  else
  {
    *p_forward = p_rates->dForward;
    *p_reverse = p_rates->dReverse;
  }

}

/*##############################################################################
// Libnucnet__Zone__getAbundances()
//############################################################################*/

gsl_vector *
Libnucnet__Zone__getAbundances(
  const Libnucnet__Zone *self
)
{

  typedef struct {
    gsl_vector *pVector;
    const Libnucnet__Zone *pZone;
  } user_data;

  gsl_vector *p_vector;

  user_data *p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !self ) LIBNUCNET__ERROR( "Invalid input" );
    
  p_user_data->pVector =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies( self->pNet->pNuc )
    );

  p_user_data->pZone = self;

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__getAbundancesCallback,
    p_user_data
  );

  p_vector = p_user_data->pVector;
  free( p_user_data );

  return p_vector;

}

/*##############################################################################
// Libnucnet__Zone__getAbundancesCallback()
//############################################################################*/

void
Libnucnet__Zone__getAbundancesCallback(
  Libnucnet__Species *p_species, void *p_data, xmlChar *sx_species
)
{

  typedef struct {
    gsl_vector *pVector;
    Libnucnet__Zone *pZone;
  } user_data;

  size_t i_index;
  user_data *p_user_data = ( user_data * ) p_data;

  if( !sx_species )
    LIBNUCNET__ERROR( "No such species" );

  i_index = Libnucnet__Species__getIndex( p_species );

  gsl_vector_set(
    p_user_data->pVector,
    i_index,
    Libnucnet__Zone__getSpeciesAbundance( p_user_data->pZone, p_species )
  );

}

/*##############################################################################
// Libnucnet__Zone__getAbundanceChanges()
//############################################################################*/

gsl_vector *
Libnucnet__Zone__getAbundanceChanges(
  const Libnucnet__Zone *self
)
{

  typedef struct {
    gsl_vector *pVector;
    const Libnucnet__Zone *pZone;
  } user_data;

  gsl_vector *p_vector;

  user_data *p_user_data = ( user_data * ) malloc( sizeof( user_data ) );
    
  if( !self ) LIBNUCNET__ERROR( "Invalid input" );
    
  p_user_data->pVector =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies( self->pNet->pNuc )
    );

  p_user_data->pZone = self;

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__getAbundanceChangesCallback,
    p_user_data
  );

  p_vector = p_user_data->pVector;
  free( p_user_data );

  return p_vector;

}

/*##############################################################################
// Libnucnet__Zone__getAbundanceChangesCallback()
//############################################################################*/

void
Libnucnet__Zone__getAbundanceChangesCallback(
  Libnucnet__Species *p_species, void *p_data, xmlChar *sx_species
)
{

  typedef struct {
    gsl_vector *pVector;
    Libnucnet__Zone *pZone;
  } user_data;

  size_t i_index;
  user_data *p_user_data = ( user_data * ) p_data;

  if( !sx_species )
    LIBNUCNET__ERROR( "No such species" );

  i_index = Libnucnet__Species__getIndex( p_species );

  gsl_vector_set(
    p_user_data->pVector,
    i_index,
    Libnucnet__Zone__getSpeciesAbundanceChange( p_user_data->pZone, p_species )
  );

}

/*##############################################################################
// Libnucnet__Zone__updateAbundances()
//############################################################################*/

void
Libnucnet__Zone__updateAbundances(
  Libnucnet__Zone *self, gsl_vector *p_abundances
)
{

  typedef struct {
    gsl_vector *pVector;
    Libnucnet__Zone *pZone;
  } user_data;

  user_data *p_user_data = ( user_data * ) malloc( sizeof( user_data ) );
    
  if(
     p_abundances->size !=
     Libnucnet__Nuc__getNumberOfSpecies(
       Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( self ) )
     )
  )
    LIBNUCNET__ERROR(
      "Number of species in vector does not match number in collection"
    );

  p_user_data->pVector = p_abundances;
  p_user_data->pZone = self;

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__updateAbundancesCallback,
    p_user_data
  ); 

  free( p_user_data );

}

/*##############################################################################
// Libnucnet__Zone__updateAbundancesCallback()
//############################################################################*/

void
Libnucnet__Zone__updateAbundancesCallback(
  Libnucnet__Species *p_species, void *p_data, xmlChar *sx_species
)
{

  typedef struct {
    gsl_vector *pVector;
    Libnucnet__Zone *pZone;
  } user_data;

  size_t i;
  user_data *p_user_data = ( user_data * ) p_data;

  if( !sx_species )
    LIBNUCNET__ERROR( "No such species" );

  i = Libnucnet__Species__getIndex( p_species );
  Libnucnet__Zone__updateSpeciesAbundance(
    p_user_data->pZone,
    p_species,
    gsl_vector_get( p_user_data->pVector, i )
  );

}
/*##############################################################################
// Libnucnet__Zone__updateAbundanceChanges()
//############################################################################*/

void
Libnucnet__Zone__updateAbundanceChanges(
  Libnucnet__Zone *self, gsl_vector *p_abundance_changes
)
{

  typedef struct {
    gsl_vector *pVector;
    Libnucnet__Zone *pZone;
  } user_data;

  user_data *p_user_data = ( user_data * ) malloc( sizeof( user_data ) );
    
  if(
     p_abundance_changes->size !=
     Libnucnet__Nuc__getNumberOfSpecies(
       Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( self ) )
     )
  )
    LIBNUCNET__ERROR(
      "Number of species in vector does not match number in collection"
    );

  p_user_data->pVector = p_abundance_changes;
  p_user_data->pZone = self;

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__updateAbundanceChangesCallback,
    p_user_data
  ); 

  free( p_user_data );

}

/*##############################################################################
// Libnucnet__Zone__updateAbundanceChangesCallback()
//############################################################################*/

void
Libnucnet__Zone__updateAbundanceChangesCallback(
  Libnucnet__Species *p_species, void *p_data, xmlChar *sx_species
)
{

  typedef struct {
    gsl_vector *pVector;
    Libnucnet__Zone *pZone;
  } user_data;

  size_t i;
  user_data *p_user_data = ( user_data * ) p_data;

  if( !sx_species )
    LIBNUCNET__ERROR( "No such species" );

  i = Libnucnet__Species__getIndex( p_species );
  Libnucnet__Zone__updateSpeciesAbundanceChange(
    p_user_data->pZone,
    p_species,
    gsl_vector_get( p_user_data->pVector, i )
  );

}

/*##############################################################################
// Libnucnet__iterateZones()
//############################################################################*/

void
Libnucnet__iterateZones(
  const Libnucnet *self,
  Libnucnet__Zone__iterateFunction pf_iterator,
  void *p_user_data
)
{

  Libnucnet__Zone **p_zones;
  size_t i, i_zones;

  if( !self ) LIBNUCNET__ERROR( "Invalid input" );

  p_zones =
    Libnucnet__createZoneArray( self );

  i_zones = Libnucnet__getNumberOfZones( self );

  for( i = 0; i < i_zones; i++ )
    if( pf_iterator( p_zones[i], p_user_data ) == 0 ) break;

  free( p_zones );
  

}

/*##############################################################################
// Libnucnet__createZoneArray().
//############################################################################*/

Libnucnet__Zone **
Libnucnet__createZoneArray(
  const Libnucnet *self
)
{

  struct {
    size_t iIndex;
    Libnucnet__Zone **pZones;
  } work_data;

  if( !self )
    LIBNUCNET__ERROR( "Invalid input" );

  work_data.pZones =
    ( Libnucnet__Zone ** )
    malloc(
      sizeof( Libnucnet__Zone ) *
      Libnucnet__getNumberOfZones( self )
    );

  if( !work_data.pZones )
    LIBNUCNET__ERROR( "Couldn't allocate memory" );

  work_data.iIndex = 0;

  xmlHashScan(
    self->pZoneHash,
    (xmlHashScanner) Libnucnet__createZoneArrayCallback,
    &work_data
  );

  if( self->pfZoneCompare )
    Libnucnet__sortZoneArray( self, work_data.pZones );

  return work_data.pZones;

}
  
/*##############################################################################
// Libnucnet__createZoneArrayCallback().
//############################################################################*/

void
Libnucnet__createZoneArrayCallback(
  Libnucnet__Zone *p_zone,
  void *p_data,
  const xmlChar *sx_zone
)
{

  typedef struct {
    size_t iIndex;
    Libnucnet__Zone **pZones;
  } work_data;

  work_data *p_work_data = ( work_data * ) p_data;

  if( !sx_zone )
    LIBNUCNET__ERROR( "No such zone" );

  p_work_data->pZones[p_work_data->iIndex++] = p_zone;

}

/*##############################################################################
// Libnucnet__sort_helper().
//############################################################################*/

int
Libnucnet__sort_helper( const void *p_1, const void *p_2 )
{

  typedef struct {
    Libnucnet__Zone *pZone;
    Libnucnet__Zone__compare_function pfFunc;
  } sort_data;

  sort_data *p_sort_1;
  sort_data *p_sort_2;

  p_sort_1 = *( sort_data * const * ) p_1;
  p_sort_2 = *( sort_data * const * ) p_2;

  return
    p_sort_1->pfFunc(
      p_sort_1->pZone,
      p_sort_2->pZone
    );

}

/*##############################################################################
// Libnucnet__sortZoneArray().
//############################################################################*/

void
Libnucnet__sortZoneArray(
  const Libnucnet *self,
  Libnucnet__Zone **p_zones
)
{

  typedef struct {
    Libnucnet__Zone *pZone;
    Libnucnet__Zone__compare_function pfFunc;
  } sort_data;

  size_t i;
  sort_data **p_sort_array;

  p_sort_array =
    (sort_data **)
    malloc(
      sizeof( sort_data * ) * Libnucnet__getNumberOfZones( self )
    );

  for( i = 0; i < Libnucnet__getNumberOfZones( self ); i++ )
  {
    p_sort_array[i] = (sort_data *) malloc( sizeof( sort_data ) );
    p_sort_array[i]->pZone = p_zones[i];
    p_sort_array[i]->pfFunc = self->pfZoneCompare;
  }

  qsort(
    p_sort_array,
    Libnucnet__getNumberOfZones( self ),
    sizeof( sort_data * ),
    Libnucnet__sort_helper
  );

  for( i = 0; i < Libnucnet__getNumberOfZones( self ); i++ )
  {
    p_zones[i] = p_sort_array[i]->pZone;
    free( p_sort_array[i] );
  }

  free( p_sort_array );

}

/*##############################################################################
// Libnucnet__writeZoneDataToXmlFile()
//############################################################################*/

void
Libnucnet__writeZoneDataToXmlFile(
  const Libnucnet *self, const char *s_output_xml_filename
)
{

  xmlDocPtr p_doc;
  xmlChar *sx_str, *sx_str2;

  /*============================================================================
  // Check input structure.
  //==========================================================================*/

  if( !self ) {
     LIBNUCNET__ERROR( "Invalid input structure" );
  }

  /*============================================================================
  // Create the xml document.
  //==========================================================================*/

  p_doc = Libnucnet__makeZoneDataXmlDocument( self );

  /*============================================================================
  // Create namespace attributes and add.
  //==========================================================================*/

  sx_str = xmlCharStrdup( W3C__NAMESPACE );
  sx_str2 = xmlCharStrdup( XSI );
  xmlNewNs( xmlDocGetRootElement( p_doc ), sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( XSI_LOCATION );
  sx_str2 = xmlCharStrdup( LIBNUCNET__ZONE_DATA__SCHEMALOCATION );
  xmlNewProp( xmlDocGetRootElement( p_doc ), sx_str, sx_str2 );
  xmlFree( sx_str2 );
  xmlFree( sx_str );

  /*============================================================================
  // Write file.
  //==========================================================================*/

  if(
       xmlSaveFormatFileEnc( s_output_xml_filename, p_doc, "UTF-8", 1 ) == -1
  ) {
      LIBNUCNET__ERROR(
        "xmlSaveFormatFileEnc: failed\n\tException #%d\n"
      );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// Libnucnet__makeZoneDataXmlDocument().
//############################################################################*/

xmlDocPtr
Libnucnet__makeZoneDataXmlDocument(
  const Libnucnet *self
)
{

  typedef struct {
    xmlNodePtr pRoot;
    xmlChar *sxMassFractionFormat;
  } work;

  work *p_work;
  xmlDocPtr p_doc;

  /*============================================================================
  // Create work structure.
  //==========================================================================*/

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCNET__ERROR( "Couldn't allocate work memory" );

  /*============================================================================
  // Get document, root, and format.
  //==========================================================================*/

  p_doc = Libnucnet__newZoneDataXmlDocument( );

  p_work->pRoot = xmlDocGetRootElement( p_doc );

  p_work->sxMassFractionFormat = self->sxMassFractionFormat;

  /*============================================================================
  // Add children.
  //==========================================================================*/

  Libnucnet__iterateZones(
    self,
    (Libnucnet__Zone__iterateFunction)
       Libnucnet__Zone__addToZoneDataXmlIterator,
    p_work
  );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  free( p_work );

  /*============================================================================
  // Return.
  //==========================================================================*/

  return p_doc;

}

/*##############################################################################
// Libnucnet__newZoneDataXmlDocument().
//############################################################################*/

xmlDocPtr
Libnucnet__newZoneDataXmlDocument( void )
{

  xmlDocPtr p_doc;
  xmlNodePtr p_root;
  xmlChar *sx_str;

  /*============================================================================
  // Get document.
  //==========================================================================*/

  sx_str = xmlCharStrdup( XML_VERSION );
  p_doc = xmlNewDoc( sx_str );
  xmlFree( sx_str );

  if (p_doc == NULL) {
    LIBNUCNET__ERROR( "DOMImplementation.createDocument: failed" );
  }

  /*============================================================================
  // Add root.
  //==========================================================================*/

  sx_str = xmlCharStrdup( ZONE_DATA );
  p_root = xmlNewNode( NULL, sx_str );
  xmlFree( sx_str );

  xmlDocSetRootElement( p_doc, p_root );

  return p_doc;

}


/*##############################################################################
// Libnucnet__Zone__addToZoneDataXmlIterator().
//############################################################################*/

int
Libnucnet__Zone__addToZoneDataXmlIterator(
  Libnucnet__Zone *self,
  void *p_data
)
{

  typedef struct {
    xmlNodePtr pRoot;
    xmlChar *sxMassFractionFormat;
  } work;

  typedef struct {
    Libnucnet__Zone *pZone;
    xmlNodePtr pOptionalPropertiesNode;
  } user_data;

  size_t i;
  Libnucnet__Species **p_species_array;
  xmlChar *sx_str;
  xmlNodePtr p_zone_node, p_mass_fractions_node;

  work *p_work = ( work * ) p_data;

  user_data *p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data ) LIBNUCNET__ERROR( "Couldn't allocate user data" );

  p_user_data->pZone = self;

  sx_str = xmlCharStrdup( ZONE );
  p_zone_node = xmlNewChild( p_work->pRoot, NULL, sx_str, NULL );
  xmlFree( sx_str );

  if( strcmp( (char *) self->sxLabel[0], ZERO ) )
  {
    sx_str = xmlCharStrdup( LABEL_1 );
    xmlNewProp( p_zone_node, sx_str, self->sxLabel[0] );
    xmlFree( sx_str );
  }

  if( strcmp( (char *) self->sxLabel[1], ZERO ) )
  {
    sx_str = xmlCharStrdup( LABEL_2 );
    xmlNewProp( p_zone_node, sx_str, self->sxLabel[1] );
    xmlFree( sx_str );
  }

  if( strcmp( (char *) self->sxLabel[2], ZERO ) )
  {
    sx_str = xmlCharStrdup( LABEL_3 );
    xmlNewProp( p_zone_node, sx_str, self->sxLabel[2] );
    xmlFree( sx_str );
  }

  /*============================================================================
  // Add optional properties.
  //==========================================================================*/

  if( xmlHashSize( self->pPropertyHash ) > 0 )
  {
    sx_str = xmlCharStrdup( OPTIONAL_PROPERTIES );
    p_user_data->pOptionalPropertiesNode =
      xmlNewChild( p_zone_node, NULL, sx_str, NULL );
    xmlFree( sx_str );

    Libnucnet__Zone__iterateOptionalProperties(
      self,
      NULL,
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function)
         Libnucnet__Zone__add_property_to_xml,
      p_user_data
    );
  }

  /*============================================================================
  // Add mass fractions.
  //==========================================================================*/

  sx_str = xmlCharStrdup( MASS_FRACTIONS );
  p_mass_fractions_node =
    xmlNewChild( p_zone_node, NULL, sx_str, NULL );
  xmlFree( sx_str );

  p_species_array =
    Libnucnet__Nuc__createSpeciesArray(
      Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( self ) )
    );

  Libnucnet__Zone__updateProperty(
    self,
    ZONE_MASS_FRACTION_FORMAT,
    NULL,
    NULL,
    (const char *) p_work->sxMassFractionFormat
  );

  for(
    i = 0;
    i < Libnucnet__Nuc__getNumberOfSpecies(
          Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( self ) )
        );
    i++
  )
    Libnucnet__Zone__addSpeciesMassFractionToXml(
      self,
      p_mass_fractions_node,
      p_species_array[i]
    );

  Libnucnet__Zone__removeProperty(
    self,
    ZONE_MASS_FRACTION_FORMAT,
    NULL,
    NULL
  );

  free( p_species_array );

  /*============================================================================
  // Clean up and return.
  //==========================================================================*/

  free( p_user_data );

  return 1;

}
  
/*##############################################################################
// Libnucnet__Zone__addSpeciesMassFractionToXml().
//############################################################################*/

void
Libnucnet__Zone__addSpeciesMassFractionToXml(
  Libnucnet__Zone *self,
  xmlNodePtr p_mass_fractions_node,
  Libnucnet__Species *p_species
)
{

  xmlNodePtr p_nuclide;
  xmlChar *sx_str1, *sx_str2;
  double d_x = 0.;

  d_x =
    Libnucnet__Zone__getSpeciesAbundance( self, p_species ) *
    Libnucnet__Species__getA( p_species );

  if( !WnMatrix__value_is_zero( d_x ) )
  {
 
    sx_str1 = xmlCharStrdup( NUCLIDE );
    p_nuclide = xmlNewChild( p_mass_fractions_node, NULL, sx_str1, NULL );
    xmlFree( sx_str1 );

    sx_str1 = xmlCharStrdup( SPECIES_NAME );
    xmlNewProp( p_nuclide, sx_str1, p_species->sxName );
    xmlFree( sx_str1 );

    sx_str1 = xmlCharStrdup( ATOMIC_NUMBER );
    sx_str2 = (xmlChar *) malloc( sizeof( xmlChar ) * NUC_BUF_SIZE );
    xmlStrPrintf(
      sx_str2,
      NUC_BUF_SIZE,
      (const WnNetChar *) "%d",
      p_species->iZ
    );
    xmlNewChild( p_nuclide, NULL, sx_str1, sx_str2 );
    xmlFree( sx_str1 );
    xmlFree( sx_str2 );

    sx_str1 = xmlCharStrdup( MASS_NUMBER );
    sx_str2 = (xmlChar *) malloc( sizeof( xmlChar ) * NUC_BUF_SIZE );
    xmlStrPrintf(
      sx_str2,
      NUC_BUF_SIZE,
      (const WnNetChar *) "%d",
      p_species->iA
    );
    xmlNewChild( p_nuclide, NULL, sx_str1, sx_str2 );
    xmlFree( sx_str1 );
    xmlFree( sx_str2 );

    sx_str1 = xmlCharStrdup( MASS_FRACTION );
    sx_str2 = (xmlChar *) malloc( sizeof( xmlChar ) * NUC_BUF_SIZE );
    xmlStrPrintf(
      sx_str2,
      NUC_BUF_SIZE,
      (const WnNetChar *)
        Libnucnet__Zone__getProperty(
          self,
          ZONE_MASS_FRACTION_FORMAT,
          NULL,
          NULL
        ),
      d_x
    );
    xmlNewChild( p_nuclide, NULL, sx_str1, sx_str2 );
    xmlFree( sx_str1 );
    xmlFree( sx_str2 );

  }

}
  
/*##############################################################################
// Libnucnet__Zone__add_property_to_xml().
//############################################################################*/

void
Libnucnet__Zone__add_property_to_xml(
  xmlChar *sx_name,
  xmlChar *sx_tag1,
  xmlChar *sx_tag2,
  xmlChar *sx_value,
  void *p_user_data
)
{

  typedef struct {
    Libnucnet__Zone *pZone;
    xmlNodePtr pOptionalPropertiesNode;
  } user_data;

  xmlNodePtr p_node;
  xmlChar *sx_str;
  
  user_data *p_data = ( user_data * ) p_user_data;

  /*============================================================================
  // Add property.
  //==========================================================================*/

  sx_str = xmlCharStrdup( PROPERTY );

  p_node =
    xmlNewChild(
      p_data->pOptionalPropertiesNode,
      NULL,
      sx_str,
      sx_value
    );

  xmlFree( sx_str );

  /*============================================================================
  // Add property attributes.
  //==========================================================================*/

  sx_str = xmlCharStrdup( PROPERTY_NAME );
  xmlNewProp( p_node, sx_str, sx_name );
  xmlFree( sx_str );

  if( sx_tag1 )
  {
    sx_str = xmlCharStrdup( PROPERTY_TAG1 );
    xmlNewProp( p_node, sx_str, sx_tag1 );
    xmlFree( sx_str );
  }

  if( sx_tag2 )
  {
    sx_str = xmlCharStrdup( PROPERTY_TAG2 );
    xmlNewProp( p_node, sx_str, sx_tag2 );
    xmlFree( sx_str );
  }

}
  
/*##############################################################################
// Libnucnet__Zone__updateProperty().
//############################################################################*/

int
Libnucnet__Zone__updateProperty(
  Libnucnet__Zone *self,
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value
)
{

  xmlChar *sx_value;

  if( !self ) LIBNUCNET__ERROR( "Updating property of invalid zone" );

  sx_value = xmlCharStrdup( s_value );

  if( !sx_value )
    LIBNUCNET__ERROR( "Couldn't allocate strings" );

  if( s_tag2 && !s_tag1 )
    LIBNUCNET__ERROR( "Expecting tag1 if tag2 present" );

  if(
      (
        xmlHashUpdateEntry3(
          self->pPropertyHash,
          (const xmlChar *) s_name,
          (const xmlChar *) s_tag1,
          (const xmlChar *) s_tag2,
          sx_value,
          (xmlHashDeallocator) xmlFree
        ) == 0
      )
  )
    return 1;
  else
  {
    xmlFree( sx_value );
    return 0;
  }

}

/*##############################################################################
// Libnucnet__Zone__removeProperty().
//############################################################################*/

int
Libnucnet__Zone__removeProperty(
  Libnucnet__Zone *self,
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2
)
{

  if( !self ) LIBNUCNET__ERROR( "Removing property from invalid zone" );

  return
    xmlHashRemoveEntry3(
      self->pPropertyHash,
      (const xmlChar *) s_name,
      (const xmlChar *) s_tag1,
      (const xmlChar *) s_tag2,
      (xmlHashDeallocator) xmlFree
    )
    + 1;

}

/*##############################################################################
// Libnucnet__Zone__getProperty().
//############################################################################*/

const char *
Libnucnet__Zone__getProperty(
  const Libnucnet__Zone *self,
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2
)
{

  if( !self ) LIBNUCNET__ERROR( "Seeking property from invalid zone" );

  if( !self->pPropertyHash )
    LIBNUCNET__ERROR( "Zone does not have a property hash" );

  return
    (const char *)
    xmlHashLookup3(
      self->pPropertyHash,
      (const xmlChar *) s_name,
      (const xmlChar *) s_tag1,
      (const xmlChar *) s_tag2
    );

}

/*##############################################################################
// Libnucnet__Zone__iterateOptionalProperties().
//############################################################################*/

void
Libnucnet__Zone__iterateOptionalProperties(
  const Libnucnet__Zone *self,
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  Libnucnet__Zone__optional_property_iterate_function pf_func,
  void *p_user_data
)
{

  typedef struct {
    Libnucnet__Zone__optional_property_iterate_function pfFunc;
    void *pUserData;
  } work_data;

  work_data *p_work;

  if( !self )
    LIBNUCNET__ERROR( "Trying to iterate properties of an invalid zone" );

  p_work = ( work_data * ) malloc( sizeof( work_data ) );

  p_work->pfFunc = pf_func;
  p_work->pUserData = p_user_data;

  xmlHashScanFull3(
    self->pPropertyHash,
    (const xmlChar *) s_name,
    (const xmlChar *) s_tag1,
    (const xmlChar *) s_tag2,
    (xmlHashScannerFull)
       Libnucnet__Zone__iterateOptionalPropertiesHelper,
    p_work
  );

  free( p_work );

}

/*##############################################################################
// Libnucnet__Zone__iterateOptionalPropertiesHelper().
//############################################################################*/

void
Libnucnet__Zone__iterateOptionalPropertiesHelper(
  const xmlChar *sx_value,
  void *p_data,
  const xmlChar *sx_name,
  const xmlChar *sx_tag1,
  const xmlChar *sx_tag2
)
{

  typedef struct {
    Libnucnet__Zone__optional_property_iterate_function pfFunc;
    void *pUserData;
  } work_data;

  work_data *p_work = (work_data *) p_data;

  p_work->pfFunc(
    (const char *) sx_name,
    (const char *) sx_tag1,
    (const char *) sx_tag2,
    (const char *) sx_value,
    p_work->pUserData
  );

}

/*##############################################################################
// Libnucnet__Zone__iterateNetViews().
//############################################################################*/

void
Libnucnet__Zone__iterateNetViews(
  const Libnucnet__Zone *self,
  const char *s_label1,
  const char *s_label2,
  const char *s_label3,
  Libnucnet__NetView__iterateFunction pf_func,
  void *p_user_data
)
{

  typedef struct {
    Libnucnet__NetView__iterateFunction pfFunc;
    void *pUserData;
  } work_data;

  work_data *p_work;

  if( !self )
    LIBNUCNET__ERROR( "Trying to iterate net views of an invalid zone" );

  p_work = ( work_data * ) malloc( sizeof( work_data ) );

  p_work->pfFunc = pf_func;
  p_work->pUserData = p_user_data;

  xmlHashScanFull3(
    self->pViewHash,
    (const xmlChar *) s_label1,
    (const xmlChar *) s_label2,
    (const xmlChar *) s_label3,
    (xmlHashScannerFull) Libnucnet__Zone__iterateNetViewsHelper,
    p_work
  );

  free( p_work );

}

/*##############################################################################
// Libnucnet__Zone__iterateNetViewsHelper().
//############################################################################*/

void
Libnucnet__Zone__iterateNetViewsHelper(
  Libnucnet__NetView * p_view,
  void *p_data,
  const xmlChar *sx_label1,
  const xmlChar *sx_label2,
  const xmlChar *sx_label3
)
{

  typedef struct {
    Libnucnet__NetView__iterateFunction pfFunc;
    void *pUserData;
  } work_data;

  work_data *p_work = (work_data *) p_data;

  p_work->pfFunc(
    p_view,
    (const char *) sx_label1,
    (const char *) sx_label2,
    (const char *) sx_label3,
    p_work->pUserData
  );

}

/*##############################################################################
// Libnucnet__relabelZone().
//############################################################################*/

int
Libnucnet__relabelZone(
  Libnucnet *self,
  Libnucnet__Zone *p_zone,
  const char *s_new_label_1,
  const char *s_new_label_2,
  const char *s_new_label_3
)
{

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !self || !p_zone )
    LIBNUCNET__ERROR( "Invalid input" );

  /*============================================================================
  // Check that zone with new labels doesn't already exist.
  //==========================================================================*/

  if(
    xmlHashLookup3(
      self->pZoneHash,
      (const xmlChar *) s_new_label_1,
      (const xmlChar *) s_new_label_2,
      (const xmlChar *) s_new_label_3
    )
  )
    return 0;

  /*============================================================================
  // Remove zone from hash.
  //==========================================================================*/

  xmlHashRemoveEntry3(
    self->pZoneHash,
    p_zone->sxLabel[0],
    p_zone->sxLabel[1],
    p_zone->sxLabel[2],
    NULL
  );

  /*============================================================================
  // Relabel.
  //==========================================================================*/

  xmlFree( p_zone->sxLabel[0] );
  xmlFree( p_zone->sxLabel[1] );
  xmlFree( p_zone->sxLabel[2] );

  if( s_new_label_1 )
    p_zone->sxLabel[0] = xmlCharStrdup( s_new_label_1 );
  else
    p_zone->sxLabel[0] = xmlCharStrdup( ZERO );

  if( s_new_label_2 )
    p_zone->sxLabel[1] = xmlCharStrdup( s_new_label_2 );
  else
    p_zone->sxLabel[1] = xmlCharStrdup( ZERO );

  if( s_new_label_3 )
    p_zone->sxLabel[2] = xmlCharStrdup( s_new_label_3 );
  else
    p_zone->sxLabel[2] = xmlCharStrdup( ZERO );

  /*============================================================================
  // Restore zone to hash.
  //==========================================================================*/

  return
    xmlHashAddEntry3(
      self->pZoneHash,
      p_zone->sxLabel[0],
      p_zone->sxLabel[1],
      p_zone->sxLabel[2],
      p_zone
    ) + 1;

}

/*##############################################################################
// Libnucnet__setZoneCompareFunction().
//############################################################################*/

void
Libnucnet__setZoneCompareFunction(
  Libnucnet *self,
  Libnucnet__Zone__compare_function pf_func
)
{

  self->pfZoneCompare = pf_func;

}

/*##############################################################################
// Libnucnet__clearZoneCompareFunction().
//############################################################################*/

void
Libnucnet__clearZoneCompareFunction(
  Libnucnet *self
)
{

  self->pfZoneCompare = NULL;

}

/*##############################################################################
// Libnucnet__Net__makeXmlDocument().
//############################################################################*/

xmlDocPtr
Libnucnet__Net__makeXmlDocument(
  const Libnucnet__Net *self
)
{

  xmlDocPtr p_doc, p_tmp_doc;
  xmlNodePtr p_root, p_node;
  xmlChar *sx_str;

  sx_str = xmlCharStrdup( XML_VERSION );
  p_doc = xmlNewDoc( sx_str );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( NUCLEAR_NETWORK );
  p_root = xmlNewNode( NULL, sx_str );
  xmlFree( sx_str );

  xmlDocSetRootElement( p_doc, p_root );

  /*============================================================================
  // Nuc.
  //==========================================================================*/

  p_tmp_doc = 
    Libnucnet__Nuc__makeXmlDocument(
      Libnucnet__Net__getNuc( self )
    );

  p_node = xmlDocGetRootElement( p_tmp_doc );

  xmlUnlinkNode( p_node );

  xmlFreeDoc( p_tmp_doc );

  xmlAddChild( p_root, p_node );

  /*============================================================================
  // Reac.
  //==========================================================================*/

  p_tmp_doc = 
    Libnucnet__Reac__makeXmlDocument(
      Libnucnet__Net__getReac( self )
    );

  p_node = xmlDocGetRootElement( p_tmp_doc );

  xmlUnlinkNode( p_node );

  xmlFreeDoc( p_tmp_doc );
  xmlCleanupParser();

  xmlAddChild( p_root, p_node );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return p_doc;

}
 
/*##############################################################################
// Libnucnet__Net__writeToXmlFile().
//############################################################################*/

void
Libnucnet__Net__writeToXmlFile(
  const Libnucnet__Net *self,
  const char *s_output_xml_filename
)
{

  xmlDocPtr p_doc;

  /*============================================================================
  // Create document.
  //==========================================================================*/

  p_doc = Libnucnet__Net__makeXmlDocument( self );

  /*============================================================================
  // Write file.
  //==========================================================================*/

  if(
       xmlSaveFormatFileEnc( s_output_xml_filename, p_doc, "UTF-8", 1 ) == -1
  ) {
      LIBNUCNET__ERROR(
        "DOMImplementation.saveDocToFile: failed\n"
      );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// Libnucnet__makeXmlDocument().
//############################################################################*/

xmlDocPtr
Libnucnet__makeXmlDocument(
  const Libnucnet *self
)
{

  xmlDocPtr p_doc, p_tmp_doc;
  xmlNodePtr p_root, p_node;
  xmlChar *sx_str;

  sx_str = xmlCharStrdup( XML_VERSION );
  p_doc = xmlNewDoc( sx_str );
  xmlFree( sx_str );

  sx_str = xmlCharStrdup( LIBNUCNET_INPUT );
  p_root = xmlNewNode( NULL, sx_str );
  xmlFree( sx_str );

  xmlDocSetRootElement( p_doc, p_root );

  /*============================================================================
  // Net.
  //==========================================================================*/

  p_tmp_doc = Libnucnet__Net__makeXmlDocument( self->pNet );

  p_node = xmlDocGetRootElement( p_tmp_doc );

  xmlUnlinkNode( p_node );

  xmlFreeDoc( p_tmp_doc );
  xmlCleanupParser();

  xmlAddChild( p_root, p_node );

  /*============================================================================
  // Zones.
  //==========================================================================*/

  p_tmp_doc = Libnucnet__makeZoneDataXmlDocument( self );

  p_node = xmlDocGetRootElement( p_tmp_doc );

  xmlUnlinkNode( p_node );

  xmlFreeDoc( p_tmp_doc );
  xmlCleanupParser();

  xmlAddChild( p_root, p_node );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return p_doc;

}
 
/*##############################################################################
// Libnucnet__writeToXmlFile().
//############################################################################*/

void
Libnucnet__writeToXmlFile(
  const Libnucnet *self,
  const char *s_output_xml_filename
)
{

  xmlDocPtr p_doc;

  /*============================================================================
  // Create document.
  //==========================================================================*/

  p_doc = Libnucnet__makeXmlDocument( self );

  /*============================================================================
  // Write file.
  //==========================================================================*/

  if(
       xmlSaveFormatFileEnc( s_output_xml_filename, p_doc, "UTF-8", 1 ) == -1
  ) {
      LIBNUCNET__ERROR(
        "DOMImplementation.saveDocToFile: failed\n"
      );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  xmlFreeDoc( p_doc );
  xmlCleanupParser();

}

/*##############################################################################
// Libnucnet__Zone__copy().
//############################################################################*/

Libnucnet__Zone *
Libnucnet__Zone__copy( const Libnucnet__Zone *self )
{

  Libnucnet__Zone *p_new_zone;
  gsl_vector *p_vector;

  if( !self ) LIBNUCNET__ERROR( "Trying to copy invalid zone" );

  p_new_zone =
    Libnucnet__Zone__new(
      Libnucnet__Zone__getNet( self ),
      Libnucnet__Zone__getLabel( self, 1 ), 
      Libnucnet__Zone__getLabel( self, 2 ), 
      Libnucnet__Zone__getLabel( self, 3 )
    );

  p_vector = Libnucnet__Zone__getAbundances( self );
  Libnucnet__Zone__updateAbundances( p_new_zone, p_vector );
  gsl_vector_free( p_vector );

  p_vector = Libnucnet__Zone__getAbundanceChanges( self );
  Libnucnet__Zone__updateAbundanceChanges( p_new_zone, p_vector );
  gsl_vector_free( p_vector );

  /*============================================================================
  // Copy properties.  Iterate rather than use hash copier since hash already
  // exists.
  //==========================================================================*/

  Libnucnet__Zone__iterateOptionalProperties(
    self,
    NULL,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       Libnucnet__Zone__property_copier,
    p_new_zone
  );

  /*============================================================================
  // Copy net views.
  //==========================================================================*/

  Libnucnet__Zone__copy_net_views( p_new_zone, self );

  return p_new_zone;

}

/*##############################################################################
// Libnucnet__Zone__property_copier().
//############################################################################*/

void
Libnucnet__Zone__property_copier(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_value,
  Libnucnet__Zone *p_zone
)
{

  Libnucnet__Zone__updateProperty(
    p_zone,
    s_name,
    s_tag1,
    s_tag2,
    s_value
  );

}

/*##############################################################################
// Libnucnet__Zone__view_copier().
//############################################################################*/

void
Libnucnet__Zone__view_copier(
  Libnucnet__NetView * p_view,
  Libnucnet__Zone * p_zone,
  xmlChar * sx_label1,
  xmlChar * sx_label2,
  xmlChar * sx_label3
)
{

  Libnucnet__Zone__updateNetView(
    p_zone,
    (const char *) sx_label1,
    (const char *) sx_label2,
    (const char *) sx_label3,
    Libnucnet__NetView__copy( p_view )
  );

}

/*##############################################################################
// Libnucnet__Zone__updateDataForUserRateFunction().
//############################################################################*/

int
Libnucnet__Zone__updateDataForUserRateFunction(
  Libnucnet__Zone * self,
  const char * s_function_key,
  void * p_extra_data
)
{

  Libnucnet__Reac__FD *p_fd;

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  p_fd =
    (Libnucnet__Reac__FD *)
    xmlHashLookup(
      self->pNet->pReac->pFDHash,
      (const xmlChar *) s_function_key
    );

  if( !p_fd )
  {
    fprintf( stderr, "Key: %s\n", s_function_key );
    LIBNUCNET__ERROR( "Function not registered" );
  }

  return
    xmlHashUpdateEntry(
      self->pUserDataHash,
      (const xmlChar *) s_function_key,
      p_extra_data,
      (xmlHashDeallocator) p_fd->pfDeallocator
    ) + 1;

}

/*##############################################################################
// Libnucnet__Zone__getDataForUserRateFunction().
//############################################################################*/

void *
Libnucnet__Zone__getDataForUserRateFunction(
  Libnucnet__Zone * self,
  const char * s_function_key
)
{

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  return
    xmlHashLookup(
      self->pUserDataHash,
      (const xmlChar *) s_function_key
    );

}

/*##############################################################################
// Libnucnet__updateZoneXmlMassFractionFormat().
//############################################################################*/

void
Libnucnet__updateZoneXmlMassFractionFormat(
  Libnucnet *self,
  const char *s_format
)
{

  if( !self ) LIBNUCNET__ERROR( "Invalid input" );

  xmlFree( self->sxMassFractionFormat );
  self->sxMassFractionFormat = xmlCharStrdup( s_format );

}

/*##############################################################################
// Libnucnet__Zone__getSummedAbundances().
//############################################################################*/

gsl_vector *
Libnucnet__Zone__getSummedAbundances(
  Libnucnet__Zone * self,
  const char * s_nucleon_type
)
{

  typedef struct {
    Libnucnet__Zone * pZone;
    xmlChar * sxNucleonType;
    gsl_vector * pResult;
  } work;

  gsl_vector * p_result;
  work * p_work;

  if( !self ) LIBNUCNET__ERROR( "Invalid zone" );

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work ) LIBNUCNET__ERROR( "Couldn't allocate memory" );

  p_work->pResult =
    gsl_vector_calloc(
      (size_t)
      Libnucnet__Nuc__getLargestNucleonNumber(
        Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( self ) ),
        s_nucleon_type
      ) + 1
    );

  p_work->pZone = self;
  p_work->sxNucleonType = xmlCharStrdup( s_nucleon_type );

  xmlHashScan(
    self->pNet->pNuc->pSpeciesHash,
    (xmlHashScanner) Libnucnet__Zone__getSummedAbundancesHelper,
    p_work
  );

  p_result = p_work->pResult;

  xmlFree( p_work->sxNucleonType );
  free( p_work );

  return p_result;

}

/*##############################################################################
// Libnucnet__Zone__getSummedAbundancesHelper().
//############################################################################*/

void
Libnucnet__Zone__getSummedAbundancesHelper(
  Libnucnet__Species * p_species,
  void * p_data,
  xmlChar * sx_name
)
{

  typedef struct {
    Libnucnet__Zone * pZone;
    xmlChar * sxNucleonType;
    gsl_vector * pResult;
  } work;

  size_t i_index;
  work * p_work = ( work * ) p_data;

  if( !sx_name ) LIBNUCNET__ERROR( "Invalid species" );

  if( strcmp( (const char *) p_work->sxNucleonType, ATOMIC_NUMBER ) == 0 )
    i_index = Libnucnet__Species__getZ( p_species ); 
  else if( strcmp( (const char *) p_work->sxNucleonType, NEUTRON_NUMBER ) == 0 )
  {
    i_index =
      Libnucnet__Species__getA( p_species ) -
      Libnucnet__Species__getZ( p_species );
  }
  else if( strcmp( (const char *) p_work->sxNucleonType, MASS_NUMBER ) == 0 )
    i_index = Libnucnet__Species__getA( p_species );
  else
    LIBNUCNET__ERROR( "Invalid nucleon type" );

  gsl_vector_set(
    p_work->pResult,
    i_index,
    gsl_vector_get( p_work->pResult, i_index )
    +
    Libnucnet__Zone__getSpeciesAbundance(
      p_work->pZone,
      p_species
    )
  );
     
}

/*##############################################################################
// Libnucnet__Zone__getEvolutionNetView().
//############################################################################*/

Libnucnet__NetView *
Libnucnet__Zone__getEvolutionNetView(
  Libnucnet__Zone * self
)
{

  Libnucnet__NetView * p_view, * p_new_view;

  p_view =
    Libnucnet__Zone__getNetView(
      self,
      EVOLUTION_NETWORK,
      NULL,
      NULL
    );

  if( p_view && !Libnucnet__NetView__wasNetUpdated( p_view ) ) return p_view;

  if( !p_view )
  {

    p_new_view =
      Libnucnet__NetView__new(
        Libnucnet__Zone__getNet( self ),
        NULL,
        NULL
      );

    Libnucnet__Zone__updateNetView(
      self,
      EVOLUTION_NETWORK,
      NULL,
      NULL,
      p_new_view
    );

    return p_new_view;

  }

  LIBNUCNET__ERROR(
    "Parent network of network view has been updated but view has not"
  );

}

/*##############################################################################
// Libnucnet__Net__populate_reaction_species_array().
//############################################################################*/

int
Libnucnet__Net__populate_reaction_species_array(
  Libnucnet__Reaction__Element * p_element,
  void * p_data
)
{

  typedef struct
  {
    Libnucnet__Species **pSpecies;
    Libnucnet__Nuc *pNuc;
    int iCount;
  } extra_data;

  extra_data * p_extra_data = ( extra_data * ) p_data;

  p_extra_data->pSpecies[p_extra_data->iCount++] =
    Libnucnet__Nuc__getSpeciesByName(
      p_extra_data->pNuc,
      Libnucnet__Reaction__Element__getName( p_element )
    );

  return 1;

}
    
/*##############################################################################
// Libnucnet__Net__createValidReactionLatexString().
//############################################################################*/

char *
Libnucnet__Net__createValidReactionLatexString(
  const Libnucnet__Net * self,
  const Libnucnet__Reaction * p_reaction
)
{

  typedef struct
  {
    Libnucnet__Nuc * pNuc;
    xmlChar * sxString;
  } work;

  work * p_work;
  xmlChar * sx_tmp;
  char * s_string;

  if( !self || !p_reaction ) LIBNUCNET__ERROR( "Invalid input" );

  if( !Libnucnet__Net__isValidReaction( self, p_reaction ) )
    return NULL;

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work )
    LIBNUCNET__ERROR( "Couldn't allocate memory for work structure" );

  p_work->pNuc = Libnucnet__Net__getNuc( self );
  p_work->sxString = NULL;

  Libnucnet__Reaction__iterateReactants(
    p_reaction,
    Libnucnet__Net__createValidReactionLatexStringHelper,
    p_work
  );

  sx_tmp = xmlCharStrdup( " \\to " );
  p_work->sxString = xmlStrcat( p_work->sxString, sx_tmp );
  xmlFree( sx_tmp );

  Libnucnet__Reaction__iterateProducts(
    p_reaction,
    Libnucnet__Net__createValidReactionLatexStringHelper,
    p_work
  );

  s_string =
    (char *)
    malloc( sizeof( char ) * ( (size_t) xmlStrlen( p_work->sxString ) + 1 ) );

  strcpy( s_string, (const char *) p_work->sxString );

  xmlFree( p_work->sxString );
  free( p_work );

  return s_string;

}
  
/*##############################################################################
// Libnucnet__Net__createValidReactionLatexStringHelper().
//############################################################################*/

int
Libnucnet__Net__createValidReactionLatexStringHelper(
  Libnucnet__Reaction__Element *p_element,
  void * p_data
)
{

  typedef struct
  {
    Libnucnet__Nuc * pNuc;
    xmlChar * sxString;
  } work;

  char * s_str;
  xmlChar *sx_plus, *sx_to, *sx_str, *sx_test;
  Libnucnet__Species * p_species;

  work * p_work = ( work * ) p_data;
 
  p_species =
    Libnucnet__Nuc__getSpeciesByName(
      p_work->pNuc,
      Libnucnet__Reaction__Element__getName( p_element )
    );

  if( p_species )
  {
    s_str = Libnucnet__Species__createLatexString( p_species );
    sx_str = xmlCharStrdup( s_str );
    free( s_str );
  }
  else
  {
    if( strcmp( GAMMA, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\gamma" );
    else if( strcmp( ELECTRON_ANTINEUTRINO, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\bar{\\nu}_e" );
    else if( strcmp( ELECTRON_NEUTRINO, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\nu_e" );
    else if( strcmp( MU_ANTINEUTRINO, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\bar{\\nu}_\\mu" );
    else if( strcmp( MU_NEUTRINO, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\nu_{\\mu}" );
    else if( strcmp( TAU_ANTINEUTRINO, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\bar{\\nu}_\\tau" );
    else if( strcmp( TAU_NEUTRINO, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\nu_{\\tau}" );
    else if( strcmp( ELECTRON, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "e^-" );
    else if( strcmp( POSITRON, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "e^+" );
    else if( strcmp( MU, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\mu" );
    else if( strcmp( ANTIMU, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\bar{\\mu}" );
    else if( strcmp( TAU, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\tau" );
    else if( strcmp( ANTITAU, (char *) p_element->sxName ) == 0 )
      sx_str = xmlCharStrdup( "\\bar{\\tau}" );
    else
      LIBNUCNET__ERROR( "Invalid other reaction element" );
  }

  sx_to = xmlCharStrdup( " \\to " );

  if( !p_work->sxString )
    p_work->sxString = xmlStrdup( sx_str );
  else {
    sx_test =
      xmlStrsub(
        p_work->sxString,
        xmlStrlen( p_work->sxString ) - xmlStrlen( sx_to ),
        xmlStrlen( sx_to )
      );
    if( !sx_test || xmlStrcmp( sx_test, sx_to ) ) {
      sx_plus = xmlCharStrdup( " + " );
      p_work->sxString = xmlStrcat( p_work->sxString, sx_plus );
      xmlFree( sx_plus );
    }
    if( sx_test ) xmlFree( sx_test );
    p_work->sxString =
      xmlStrcat( p_work->sxString, sx_str );
  }

  xmlFree( sx_to );
  xmlFree( sx_str );

  return 1;

}

/*##############################################################################
// Libnucnet__Zone__toggleReverseRateDetailedBalance().
//############################################################################*/

void
Libnucnet__Zone__toggleReverseRateDetailedBalance(
  Libnucnet__Zone * self,
  const char * s_switch
)
{

  if( strcmp( s_switch, LIBNUCNET_ON ) == 0 )
    self->iComputeReverseDetailedBalance = 1;
  else if( strcmp( s_switch, LIBNUCNET_OFF ) == 0 )
    self->iComputeReverseDetailedBalance = 0;
  else
    LIBNUCNET__ERROR( "Invalid switch for detailed balance" );

}

/*##############################################################################
// Libnucnet__Zone__copy_net_views().
//############################################################################*/

void
Libnucnet__Zone__copy_net_views(
  Libnucnet__Zone * p_new_zone,
  const Libnucnet__Zone * p_old_zone
)
{

  if( !p_old_zone || !p_new_zone )
    LIBNUCNET__ERROR( "Invalid zone" );

  if( p_new_zone->pNet != p_old_zone->pNet )
    LIBNUCNET__ERROR( "Underlying networks of the two zones are not the same" );

  xmlHashScanFull(
    p_old_zone->pViewHash,
    (xmlHashScannerFull) Libnucnet__Zone__view_copier,
    p_new_zone
  );

}

/*##############################################################################
// Libnucnet__Zone__clearNetViews()
//############################################################################*/

int
Libnucnet__Zone__clearNetViews( Libnucnet__Zone * self )
{

  if( !self ) return 0;

  xmlHashFree(
    self->pViewHash,
    (xmlHashDeallocator) Libnucnet__NetView__free
  );

  self->pViewHash = xmlHashCreate( 0 );

  return 1;

}

/*##############################################################################
// Libnucnet__Zone__clearRates()
//############################################################################*/

int
Libnucnet__Zone__clearRates( Libnucnet__Zone * self )
{

  if( !self ) return 0;

  xmlHashFree(
    self->pRatesHash,
    (xmlHashDeallocator) Libnucnet__Rates__free
  );

  self->pRatesHash = xmlHashCreate( 0 );

  return 1;

}

/*##############################################################################
// Libnucnet__Zone__isComputingReverseRatesFromDetailedBalance().
//############################################################################*/

int
Libnucnet__Zone__isComputingReverseRatesFromDetailedBalance(
  const Libnucnet__Zone * p_zone
)
{

  return p_zone->iComputeReverseDetailedBalance;

}
