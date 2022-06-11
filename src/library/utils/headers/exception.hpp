#pragma once

/**
 * \file exception.hpp
 *
 * \brief Definition of class \c Exception, which offers capability to
 * create exceptions  with a bit  more information on the  location of
 * the error.
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

#include <exception>
#include <iostream>
#include <string>

#define ERROR_UNLESS(expr, msg)                                                \
  if (!(expr)) {                                                               \
    throw utils::Exception(__FILE__, __LINE__, msg);                           \
  }

#define DISPLAY_EXCEPTION_INFO(e)                                              \
  cout << "Exception: " << e.get_description() << endl                         \
       << "File: " << e.get_filename() << endl                                 \
       << "Line: " << e.get_line_number() << endl;

/**
 * \defgroup UTILSNameSpace Namespace utils.
 * @{
 */

/**
 * \namespace utils
 *
 * \brief Contains definition and implementation of utility classes.
 */

namespace utils {

class Exception : public std::exception {
public:
  /**
   * \fn Exception(const std::string& file, int line, const string& desc)
   *
   * \brief Creates an instance of this class.
   *
   * \param file Name of the file  containing the source code.
   * \param line Number of the code line where the exception occurred.
   * \param desc A message describing the exception cause.
   */
  Exception(const std::string &file = "Unknown", int line = 0,
            const std::string &desc = "None");

  /** Returns a description of the cause of the exception. */
  const char *get_description() const noexcept;
  const char *what() const noexcept override;

  /**
   * Returns the name of the file containing the source code that threw the
   * exception.
   */
  const char *get_filename() const noexcept;

  /**
   * Returns the number of the source code line (in the file) that threw the
   * exception.
   */
  int get_line_number() const noexcept;

private:
  std::string
      _filename; ///< Filename of the source code that threw the exception.
  int _linenum;  ///< Number of the source code line that threw the exception.
  std::string
      _description; ///< A short description of the cause of the exception.
};

} // namespace utils
