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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief Code for wrapping zones,
//!        reactions, species, reactants and products, and properties.
//!
////////////////////////////////////////////////////////////////////////////////

#include "nnt/wrappers.hpp"

/**
 * @brief The NucNet Tools namespace.
 */
namespace nnt
{

//##############################################################################
// ReactionElement methods.
//##############################################################################

/**
 * A method that sets the wrapped Libnucnet__Reaction__Element.
 * \param p_reaction_element A pointer to the element to be wrapped.
 */

void
ReactionElement::setNucnetReactionElement(
  Libnucnet__Reaction__Element * p_reaction_element )
{
  pReactionElement = p_reaction_element;
}

/**
 * A method that returns the wrapped Libnucnet__Reaction__Element.
 * \return A pointer to the wrapped Libnucnet__Reaction__Element.
 */

Libnucnet__Reaction__Element * ReactionElement::getNucnetReactionElement()
{
  return pReactionElement;
}

//##############################################################################
// Reaction methods.
//##############################################################################

/**
 * A method that sets the wrapped Libnucnet__Reaction.
 * \param p_reaction A pointer to the Libnucnet__Reaction to be wrapped.
 */

void Reaction::setNucnetReaction( Libnucnet__Reaction * p_reaction )
{
  pReaction = p_reaction;
}

/**
 * A method that returns the wrapped Libnucnet__Reaction.
 * \return A pointer to the wrapped Libnucnet__Reaction.
 */

Libnucnet__Reaction * Reaction::getNucnetReaction()
{
  return pReaction;
}

//##############################################################################
// Species methods.
//##############################################################################

/**
 * A method that sets the wrapped Libnucnet__Species.
 * \param p_species A pointer to the Libnucnet__Species to be wrapped.
 */

void Species::setNucnetSpecies( Libnucnet__Species * p_species )
{
  pSpecies = p_species;
}

/**
 * A method that returns the wrapped Libnucnet__Species.
 * \return A pointer to the wrapped Libnucnet__Species.
 */

Libnucnet__Species * Species::getNucnetSpecies()
{
  return pSpecies;
}

//##############################################################################
// Zone methods.
//##############################################################################

/**
 * A method that sets the wrapped Libnucnet__Zone.
 * \param p_zone A pointer to the Libnucnet__Zone to be wrapped.
 */

void Zone::setNucnetZone( Libnucnet__Zone * p_zone )
{
  pZone = p_zone;
}

/**
 * A method that returns the wrapped Libnucnet__Zone.
 * \return A pointer to the wrapped Libnucnet__Zone.
 */

Libnucnet__Zone * Zone::getNucnetZone()
{
  return pZone;
}

/**
 * A method to return a Libnucnet__NetView from a zone.  If the view
 *   does not exist or if the underlying network has been updated, the
 *   routine creates, stores, and returns a new view.  To create the
 *   new view, the routine will use the first label as the nuclear
 *   XPath expression and the second as the reaction XPath expression.
 * \param s_label1 The first label for the view.
 * \param s_label2 The second label for the view (optional).
 * \param s_label3 The third label for the view (optional).
 */

Libnucnet__NetView *
Zone::getNetView(
  const char * s_label1,
  const char * s_label2,
  const char * s_label3
)
{

  Libnucnet__NetView * p_view, * p_new_view;

  p_view =
    Libnucnet__Zone__getNetView(
      this->getNucnetZone(),
      s_label1,
      s_label2,
      s_label3
    );

  if( p_view && !Libnucnet__NetView__wasNetUpdated( p_view ) ) return p_view;

  p_new_view =
    Libnucnet__NetView__new(
      Libnucnet__Zone__getNet( this->getNucnetZone() ),
      s_label1,
      s_label2
    );

  Libnucnet__Zone__updateNetView(
    this->getNucnetZone(),
    s_label1,
    s_label2,
    s_label3,
    p_new_view
  );

  return p_new_view; 

}

Libnucnet__NetView *
Zone::getNetView(
  const char * s_label1,
  const char * s_label2
)
{

  return this->getNetView( s_label1, s_label2, NULL );

}

Libnucnet__NetView *
Zone::getNetView(
  const char * s_label1
)
{

  return this->getNetView( s_label1, NULL );

}

//############################################################################
// Zone::hasProperty().
//############################################################################

/**
 * \brief Determines whether a property is present in a zone.
 *
 * \param s_name A string giving the name of the property.
 * \param s_tag1 A string giving an 
 *                 extra tag associated with the property (optional).
 * \param s_tag2 A string giving a second extra tag associated with the
 *		      property (optional).
 * \return true if the property is present and false if not.
 */

bool
Zone::hasProperty(
  const std::string s_name,
  const std::string s_tag1,
  const std::string s_tag2
)
{

  if(
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      s_tag2.c_str()
    )
  )
    return true;
  else
    return false;

} 

bool
Zone::hasProperty(
  const std::string s_name,
  const std::string s_tag1
)
{

  if(
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      NULL
    )
  )
    return true;
  else
    return false;

} 

bool
Zone::hasProperty(
  const std::string s_name
)
{

  if(
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      NULL,
      NULL
    )
  )
    return true;
  else
    return false;

} 

//############################################################################
// Zone::hasFunction().
//############################################################################

/**
 * \brief Determines whether a zone has the function with the given key.
 *
 * \param s_key A string giving the key.
 * \return true if the function is present and false if not.
 */

bool
Zone::hasFunction( std::string s_key )
{

  if(
    this->function_map.find( s_key ) != this->function_map.end()
  )
    return true;
  else
    return false;

} 

//############################################################################
// Zone::updateFunction().
//############################################################################

/**
 * \brief Update the function for a zone has the function with the given key.
 *
 * \param s_key A string giving the key.
 * \param new_function The function to update with.
 * \param s_doc A string giving documentation for the function (optional but
 *              required if s_tag present).
 * \param s_tag A string giving the tag for the function (optional).
 * \return true if the update succeeded and false if not.
 */

void
Zone::updateFunction(
  const std::string s_key,
  boost::any new_function,
  const std::string s_doc,
  const std::string s_tag
)
{

  function_map_t::iterator it = this->function_map.find( s_key );

  if( it != this->function_map.end() )
  {
    this->function_map.erase( it );
  }

  this->function_map.insert(
    zone_function(
      s_key,
      new_function,
      s_doc,
      s_tag
    )
  );

} 

void
Zone::updateFunction(
  const std::string s_key,
  boost::any new_function,
  const std::string s_doc
)
{
  this->updateFunction( s_key, new_function, s_doc, "" );
}

void
Zone::updateFunction(
  const std::string s_key,
  boost::any new_function
)
{
  this->updateFunction( s_key, new_function, "", "" );
}

//############################################################################
// Zone::getFunction().
//############################################################################

/**
 * \brief Retrieve the function with the given key.
 *
 * \param s_key A string giving the key.
 * \return The function with the given key or exit with an
 *         error if there is no function.
 */

boost::any
Zone::getFunction( const std::string s_key )
{

  function_map_t::iterator it = this->function_map.find( s_key );

  if( it != this->function_map.end() )
  {
    return it->get_func();
  }
  else
  {
    std::cerr << "No function: " << s_key << "." << std::endl;
    exit( EXIT_FAILURE );
  }  

} 

//############################################################################
// Zone::getFunctionDoc().
//############################################################################

/**
 * \brief Retrieve the documentation for the function with the given key.
 *
 * \param s_key A string giving the key.
 * \return The documentation for the function with the given key or exit with an
 *         error if there is no function.
 */

std::string
Zone::getFunctionDoc( const std::string s_key )
{

  function_map_t::iterator it = this->function_map.find( s_key );

  if( it != this->function_map.end() )
  {
    return it->get_doc();
  }
  else
  {
    std::cerr << "No function: " << s_key << "." << std::endl;
    exit( EXIT_FAILURE );
  }  

} 

//############################################################################
// Zone::getFunctionTag().
//############################################################################

/**
 * \brief Retrieve the tag for the function with the given key.
 *
 * \param s_key A string giving the key.
 * \return The tag for the function with the given key or exit with an
 *         error if there is no function.
 */

std::string
Zone::getFunctionTag( const std::string s_key )
{

  function_map_t::iterator it = this->function_map.find( s_key );

  if( it != this->function_map.end() )
  {
    return it->get_tag();
  }
  else
  {
    std::cerr << "No function: " << s_key << "." << std::endl;
    exit( EXIT_FAILURE );
  }  

} 

//############################################################################
// Zone::getListOfFunctions().
//############################################################################

/**
 * \brief Retrieves the currently defined zone function keys.
 *
 * \param s_tag A string giving the tag for the functions (optional);
 * \return A vector containing the keys of the currently defined functions;
 */

std::vector<std::string>
Zone::getListOfFunctions( const std::string s_tag )
{

  std::vector<std::string> keys;

  for(
    function_map_t::iterator it = this->function_map.begin();
    it != this->function_map.end();
    it++
  )
  {

    if( it->get_tag() == s_tag )
      keys.push_back( it->get_quantity() );

  }

  return keys;

} 

std::vector<std::string>
Zone::getListOfFunctions( )
{

  std::vector<std::string> keys;

  for(
    function_map_t::iterator it = this->function_map.begin();
    it != this->function_map.end();
    it++
  )
  {

    keys.push_back( it->get_quantity() );

  }

  return keys;

} 

} //namespace nnt
