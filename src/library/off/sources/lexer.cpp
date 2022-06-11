/**
 * \file lexer.cpp
 *
 * \brief  Implementation  of  the  class Lexer,  which  represents  a
 * lexical analyzer for scanning tokens  from an OFF file describing a
 * triangle mesh.
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

#include "lexer.hpp"

#include <exception.hpp>

#include <cctype>
#include <cstdlib>

namespace off {

Lexer::Lexer(std::istream &is) : _realstr(is) {
  _line_counter = 0;
  _char_counter = 0;
}

bool Lexer::get_string(std::string &s) {
  if (!skip_space()) {
    return false;
  }

  char c;
  if (!_realstr.get(c)) {
    return false;
  }

  if (isalpha(c)) {
    char str[256];
    int i = 0;
    do {
      str[i] = c;

      if (!_realstr.get(c)) {
        return false;
      }

      increment_char();

      ++i;
    } while (isalnum(c) && (i < 255));

    str[i] = '\0';

    ERROR_UNLESS(!isalnum(c), "Expected to read in an alphanumeric character");

    s = str;
    _realstr.putback(c);
  } else {
    _realstr.putback(c);

    return false;
  }

  return true;
}

bool Lexer::get_integer(int &x) {
  if (!skip_space()) {
    return false;
  }

  char c;
  if (!_realstr.get(c)) {
    return false;
  }

  if (isdigit(c)) {
    char str[256];
    int i = 0;
    do {
      str[i] = c;

      if (!_realstr.get(c)) {
        return false;
      }

      increment_char();

      ++i;
    } while (isdigit(c) && (i < 255));

    str[i] = '\0';

    ERROR_UNLESS(!isdigit(c), "Expected to read in a digit");

    x = atoi(str);

    _realstr.putback(c);
  } else {
    _realstr.putback(c);

    return false;
  }

  return true;
}

bool Lexer::get_double(double &d) {
  if (!skip_space()) {
    return false;
  }

  char c;
  if (!_realstr.get(c)) {
    return false;
  }

  _realstr.putback(c);

  ERROR_UNLESS((_realstr >> d), "Expected to read in a floating-point number");

  return true;
}

int Lexer::get_line_counter() const { return _line_counter; }

int Lexer::get_char_counter() const { return _char_counter; }

void Lexer::increment_line() {
  ++_line_counter;
  reset_char_counter();
}

void Lexer::increment_char() { ++_char_counter; }

void Lexer::decrement_char() {
  if (_char_counter > 0)
    --_char_counter;
}

void Lexer::reset_char_counter() { _char_counter = 0; }

bool Lexer::skip_space() {
  char c;

  for (;;) {
    if (!_realstr.get(c)) {
      return false;
    }

    increment_char();

    if (c == '\n') {
      increment_line();
    }

    if (!isspace(c)) {
      _realstr.putback(c);
      decrement_char();
      return true;
    }
  }

  return true;
}

} // namespace off
