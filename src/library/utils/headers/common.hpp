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
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#include <string>

namespace utils {

/**
 * \fn get_or_default(const M &map, const K &key, const V &value)
 *
 * \brief Returns  the value associated with  a given key in  a map if
 * the key is present, and returns a given default value otherwise.
 *
 * \param map A dictionary of pair of the form (key, value).
 * \param key A key whose associated value is to be searched for.
 * \param  value A  value to  return in  case the  key is  not in  the
 * dictionary.
 * \return The  value associated with  the key  (if any) or  the given
 * value.
 */
template <typename M, typename K, typename V>
const V &get_or_default(const M &map, const K &key, const V &value) {
  const auto it = map.find(key);
  if (it != map.end())
    return it->second;

  return value;
}

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
