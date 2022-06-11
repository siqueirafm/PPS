#pragma once

/**
 * \file common.hpp
 *
 * \brief Declaration of commonly used utility functions.
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

#include <string>

namespace utils {

/**
 * \fn to_lower(const std::string &str)
 *
 * \brief Returns a lower-case copy of a given string.
 *
 * \param str A given string.
 * \return A lower-case copy of the given string.
 */
std::string to_lower(const std::string &str);

} // namespace utils
