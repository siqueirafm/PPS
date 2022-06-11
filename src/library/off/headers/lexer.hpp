#pragma once

/**
 * \file lexer.hpp
 *
 * \brief Definition  of the class  Lexer, which represents  a lexical
 * analyzer for scanning tokens from an OFF file describing a triangle
 * mesh.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 2.0
 * \date July 2009
 */

#include <istream>
#include <string>

/**
 * \defgroup OFFNameSpace Namespace off.
 * @{
 */

/**
 * \namespace off
 *
 * \brief   The  namespace   \c  off   contains  the   definition  and
 * implementation of  classes for reading and  writing the topological
 * and geometric information of triangle meshes from / to OFF files.
 */

namespace off {

/**
 * \class Lexer
 *
 * \brief  This class  represents  a lexical  analyzer for  scanning
 * tokens from an OFF file describing a triangle mesh.
 */
class Lexer {
public:
  /**
   * \fn Lexer( std::istream& is )
   *
   * \brief Creates an instance of this class.
   *
   * \param is An input stream for an OFF file.
   */
  Lexer(std::istream &is);

  /**
   * \fn bool get_string( std::string& s )
   *
   * \brief Reads in a string from the input stream.
   *
   * \param s A reference to a string.
   */
  bool get_string(std::string &s);

  /**
   * \fn bool get_integer( int& x )
   *
   * \brief Reads in a non-negative integer from the input stream.
   *
   * \param x A reference to an integer.
   */
  bool get_integer(int &x);

  /**
   * \fn bool get_double( double& d )
   *
   * \brief  Reads  in a  double  precision  number  from the  input
   * stream.
   *
   * \param d A reference to a double precision number.
   */
  bool get_double(double &d);

  /**
   * \fn int get_line_counter() const
   *
   * \brief Returns the value of the file line counter.
   *
   * \return The value of the file line counter.
   */
  int get_line_counter() const;

  /**
   * \fn int get_char_counter() const
   *
   * \brief Returns the value of the file char counter.
   *
   * \return The value of the file char counter.
   */
  int get_char_counter() const;

private:
  /**
   * \fn void increment_line()
   *
   * \brief Increments the file line counter.
   */
  void increment_line();

  /**
   * \fn void increment_char()
   *
   * \brief Increments the file char counter.
   */
  void increment_char();

  /**
   * \fn void decrement_char()
   *
   * \brief Decrements the file char counter.
   */
  void decrement_char();

  /**
   * \fn void reset_char_counter()
   *
   * \brief Assigns the value zero to the file char counter.
   */
  void reset_char_counter();

  /**
   * \fn bool skip_space()
   *
   * \brief Reads  in the  input stream until  a character  or digit
   * shows up.
   */
  bool skip_space();

  // ---------------------------------------------------------------
  //
  // Private data members
  //
  // ---------------------------------------------------------------

  std::istream &_realstr; ///< The input stream.
  int _line_counter;      ///< The file line counter.
  int _char_counter;      ///< The file char counter.
};

} // namespace off

/** @} */ // end of group class.
