/**
 * \file exception.cpp
 *
 * \brief Implememtation of class \c Exception.
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

#include "exception.hpp"

namespace utils {
Exception::Exception(const std::string &file, int line, const std::string &desc)
    : _filename(file), _linenum(line), _description(desc) {}

const char *Exception::get_filename() const noexcept {
  return _filename.c_str();
}

int Exception::get_line_number() const noexcept { return _linenum; }

const char *Exception::get_description() const noexcept {
  return _description.c_str();
}

const char *Exception::what() const noexcept { return get_description(); }
} // namespace utils
