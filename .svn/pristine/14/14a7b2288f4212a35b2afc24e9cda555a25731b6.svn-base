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
//! \brief A header file to define classes to wrap zones,
//!        reactions, species, reactants and products, and properties.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_WRAPPERS_H
#define NNT_WRAPPERS_H

#include <iostream>
#include <string>
#include <vector>

#include <boost/bind.hpp>
#include <boost/any.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include <Libnucnet.h>
#include <Libstatmech.h>
#include <Libnuceq.h>

/**
 * @brief The NucNet Tools namespace.
 */
namespace nnt
{

  //############################################################################
  // Typedefs. 
  //############################################################################

  class zone_function
  {

    public:

      zone_function(
	const std::string& _quantity,
	boost::any _func,
        const std::string& _doc,
        const std::string& _tag
      ) :
      func( _func ), quantity( _quantity ), doc( _doc ), tag( _tag ){}

      std::string get_quantity() const { return quantity; }
      std::string get_doc() const { return doc; }
      std::string get_tag() const { return tag; }
      boost::any get_func() const { return func; }

    private:
      boost::any func;
      std::string quantity;
      std::string doc;
      std::string tag;

  };
      
  typedef boost::multi_index_container<
    zone_function,
    boost::multi_index::indexed_by<
      boost::multi_index::hashed_unique<
	boost::multi_index::const_mem_fun<
	  zone_function,
	  std::string,
	  &zone_function::get_quantity
	>
      >,
      boost::multi_index::hashed_non_unique<
	boost::multi_index::const_mem_fun<
	  zone_function,
	  std::string,
	  &zone_function::get_tag
	>
      >
    >
  > function_map_t;

  //############################################################################
  // ReactionElement.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Reaction__Element. Member functions are
   * defined in wrappers.cpp.
   */

  class ReactionElement
  {

    public:
      Libnucnet__Reaction__Element * getNucnetReactionElement();
      void setNucnetReactionElement( Libnucnet__Reaction__Element * );

    protected:
      Libnucnet__Reaction__Element * pReactionElement;

  };

  //############################################################################
  // Reaction.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Reaction.
   */

  class Reaction
  {

    public:
      Libnucnet__Reaction * getNucnetReaction();
      void setNucnetReaction( Libnucnet__Reaction * );

    protected:
      Libnucnet__Reaction * pReaction;

  };

  //############################################################################
  // Zone.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Zone.
   */

  class Zone
  {

    public:
      Zone(){}
      Libnucnet__Zone * getNucnetZone();
      void setNucnetZone( Libnucnet__Zone * );
      template<class R> R getProperty( const std::string );
      template<class R> R getProperty( const std::string, const std::string );
      template<class R> R 
         getProperty( const std::string, const std::string, const std::string );
      template<class T> int updateProperty( const std::string, T );
      template<class T>
        int updateProperty(
          const std::string, const std::string, T
        );
      template<class T>
        int updateProperty(
          const std::string, const std::string, const std::string, T
        );
      bool hasProperty( const std::string );
      bool hasProperty( const std::string, const std::string );
      bool hasProperty(
        const std::string, const std::string, const std::string
      );
      bool hasFunction( std::string );
      void updateFunction(
        const std::string,
        boost::any,
        const std::string,
        const std::string
      );
      void updateFunction( const std::string, boost::any, const std::string );
      void updateFunction( const std::string, boost::any );
      boost::any getFunction( const std::string );
      std::string getFunctionDoc( const std::string );
      std::string getFunctionTag( const std::string );
      std::vector<std::string> getListOfFunctions( const std::string );
      std::vector<std::string> getListOfFunctions( );
      Libnucnet__NetView *
        getNetView( const char *, const char *, const char * );
      Libnucnet__NetView * getNetView( const char *, const char * );
      Libnucnet__NetView * getNetView( const char * );

      bool operator==(const Zone &p) const
      {

        std::string s11 = Libnucnet__Zone__getLabel( pZone, 1 );
        std::string s12 = Libnucnet__Zone__getLabel( pZone, 2 );
        std::string s13 = Libnucnet__Zone__getLabel( pZone, 3 );

        std::string s21 = Libnucnet__Zone__getLabel( p.pZone, 1 );
        std::string s22 = Libnucnet__Zone__getLabel( p.pZone, 2 );
        std::string s23 = Libnucnet__Zone__getLabel( p.pZone, 3 );

        return (s11 == s21) && (s12 == s22) && (s13 == s23 );

      }

      bool operator<(const Zone &p) const
      {

        std::string s11 = Libnucnet__Zone__getLabel( pZone, 1 );
        std::string s12 = Libnucnet__Zone__getLabel( pZone, 2 );
        std::string s13 = Libnucnet__Zone__getLabel( pZone, 3 );

        std::string s21 = Libnucnet__Zone__getLabel( p.pZone, 1 );
        std::string s22 = Libnucnet__Zone__getLabel( p.pZone, 2 );
        std::string s23 = Libnucnet__Zone__getLabel( p.pZone, 3 );

        if( s11 != s21 )
          return s11 < s21;

        if( s12 != s22 )
          return s22 < s22;

        return s13 < s23;

      }

    protected:
      Libnucnet__Zone * pZone;
      function_map_t function_map;

  };

  //############################################################################
  // Species.
  //############################################################################

  /**
   * A class that wraps a Libnucnet__Species.
   */

  class Species
  {

    public:
      Libnucnet__Species * getNucnetSpecies();
      void setNucnetSpecies( Libnucnet__Species * );

    protected:
      Libnucnet__Species * pSpecies;

  };

//############################################################################
// Zone::updateProperty().
//############################################################################

/**
 * \brief Updates a property in a zone.
 *
 * \param s_name A string giving the name of the property.
 * \param s_tag1 A string giving an 
 *                 extra tag associated with the property (optional).
 * \param s_tag2 A string giving a second extra tag associated with the
 *		      property (optional).
 * \param value The value of the property to be updated.
 * \return 1 (true) if update succeeded or 0 (false) if not.
 */

template<class T>
int
Zone::updateProperty(
  const std::string s_name,
  const std::string s_tag1,
  const std::string s_tag2,
  T value
)
{

  return
    Libnucnet__Zone__updateProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      s_tag2.c_str(),
      boost::lexical_cast<char *>( value )
    );

} 

template <class T>
int
Zone::updateProperty(
  const std::string s_name,
  const std::string s_tag1,
  T value
)
{

  return
    Libnucnet__Zone__updateProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      s_tag1.c_str(),
      NULL,
      boost::lexical_cast<std::string>( value ).c_str()
    );

} 

template<class T>
int
Zone::updateProperty(
  const std::string s_name,
  T value
)
{

  return
    Libnucnet__Zone__updateProperty(
      this->getNucnetZone(),
      s_name.c_str(),
      NULL,
      NULL,
      boost::lexical_cast<std::string>( value ).c_str()
    );

} 

//############################################################################
// Zone::getProperty().
//############################################################################

/**
  \brief Returns the current value of a property in a zone.

  \param s_name A string giving the name of the property.
  \param s_tag1 A string giving an
		  extra tag associated with the property (optional).
  \param s_tag2 A string giving a second extra tag associated with the
		    property (optional).
  \return The value of the property.
*/

template<class R>
R
Zone::getProperty(
  const std::string s_name,
  const std::string s_tag1,
  const std::string s_tag2
)
{

  const char * s_value =
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(), s_name.c_str(), s_tag1.c_str(), s_tag2.c_str()
    );

  if( s_value )
    return boost::lexical_cast<R>( s_value );
  else
  {

    std::cerr << "Property " <<
      s_name << " " << s_tag1 <<
      " not present." << std::endl;

    exit( EXIT_FAILURE );

  }

}

template<class R>
R
Zone::getProperty(
  const std::string s_name,
  const std::string s_tag1
)
{

  const char * s_value =
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(), s_name.c_str(), s_tag1.c_str(), NULL
    );

  if( s_value )
    return boost::lexical_cast<R>( s_value );
  else
  {

    std::cerr << "Property " <<
      s_name << " " << s_tag1 <<
      " not present." << std::endl;

    exit( EXIT_FAILURE );

  }

}

template<class R>
R
Zone::getProperty(
  const std::string s_name
)
{

  const char * s_value =
    Libnucnet__Zone__getProperty(
      this->getNucnetZone(), s_name.c_str(), NULL, NULL
    );

  if( s_value )
    return boost::lexical_cast<R>( s_value );
  else
  {

    std::cerr << "Property " <<
      s_name <<
      " not present." << std::endl;

    exit( EXIT_FAILURE );

  }

}

} // namespace nnt

#endif // NNT_WRAPPERS_H
