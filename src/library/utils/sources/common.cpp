/**
 * \file common.cpp
 *
 * \brief Implementation of commonly used utility functions.
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

#include "common.hpp"

#include <algorithm>

namespace utils {

std::string to_lower(const std::string &str) {
  std::string copyStr(str);

  std::transform(std::begin(copyStr), std::end(copyStr), std::begin(copyStr),
                 [](char ch) { return std::tolower(ch); });

  return copyStr;
}

} // namespace utils
