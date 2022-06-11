/**
 * \file reader.cpp
 *
 * \brief Implementation of the  class Reader, which represents a file
 * reader for reading  in a triangle surface mesh  information from an
 * OFF file.
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

#include "reader.hpp"

#include <common.hpp>    // to_lower()
#include <exception.hpp> // ERROR_UNLESS

#include <iostream>

namespace off {

using std::cout;
using std::endl;

using utils::to_lower;

namespace fs = std::filesystem;

Reader::Reader(const fs::path &filepath) : _filepath(filepath) {
  ERROR_UNLESS(fs::exists(_filepath),
               "Given input file does not exist: " + _filepath.string());
  ERROR_UNLESS(fs::is_regular_file(_filepath),
               "Given input file path does not correspond to a file: " +
                   _filepath.string());
  const auto extension = to_lower(_filepath.extension().string());
  ERROR_UNLESS(extension == ".off",
               "Given input file name is expected to have extension '.off'");
}

void Reader::read(unsigned &nv, double *&vset, unsigned &nf, unsigned *&fset) {
  std::filebuf fb;

  const auto isOpen = fb.open(_filepath.string().c_str(), std::ios::in);
  ERROR_UNLESS(isOpen, "Could not open file " + _filepath.string());

  std::istream istr(&fb);

  Lexer lexer(istr);

  vset = 0;
  fset = 0;

  read_header(lexer, nv, nf);

  read_vertices(lexer, nv, vset);

  read_faces(lexer, nf, nv, fset);

  fb.close();
}

void Reader::read_header(Lexer &lexer, unsigned &nv, unsigned &nf) {
  // The first line should start with the string "OFF".
  std::string stroff;

  ERROR_UNLESS(lexer.get_string(stroff),
               "Error while trying to read in a string");

  ERROR_UNLESS(stroff == "OFF", "Expected to read in the word OFF");

  // Read in the number of vertices of the mesh.
  int number;

  ERROR_UNLESS(lexer.get_integer(number),
               "Error while trying to read in the number of vertices");

  nv = unsigned(number);

  ERROR_UNLESS(nv >= 4, "Number of vertices must be at least four");

  // Read in the number of faces of the mesh.
  ERROR_UNLESS(lexer.get_integer(number),
               "Error while trying to read in the number of faces");

  nf = unsigned(number);

  ERROR_UNLESS(nf >= 4, "Number of faces must be at least four");

  // Read in the number of edges of the surface.
  ERROR_UNLESS(lexer.get_integer(number),
               "Error while trying to read in the number of edges");
}

void Reader::read_vertices(Lexer &lexer, unsigned nv, double *&vset) {
  vset = new double[3 * nv];

  // Read "nv" lines with the vertex coordinates
  for (unsigned i = 0; i < nv; i++) {
    // Read in the first coordinate.
    double x;

    ERROR_UNLESS(lexer.get_double(x),
                 "Error while trying to read in vertex x-coordinate");

    // Read in the second coordinate.
    double y;

    ERROR_UNLESS(lexer.get_double(y),
                 "Error while trying to read in vertex y-coordinate");

    // Read in the third coordinate.
    double z;

    ERROR_UNLESS(lexer.get_double(z),
                 "Error while trying to read in vertex z-coordinate");

    const unsigned j = 3 * i;

    vset[j] = x;
    vset[j + 1] = y;
    vset[j + 2] = z;
  }
}

void Reader::read_faces(Lexer &lexer, unsigned nf, unsigned nv,
                        unsigned *&fset) {
  // Allocate memory for storing the vertices.
  fset = new unsigned[3 * nf];

  // Read "nf" lines of the form "3 v1 v2 v3", where "v1", "v2", and
  // "v3"  are the  vertex  indices  of the  face  defined by  these
  // vertices.
  for (unsigned i = 0; i < nf; i++) {
    // Read in the number of vertices of the i-th face.
    int number;

    ERROR_UNLESS(lexer.get_integer(number),
                 "Error while trying to read in number of face vertices");

    ERROR_UNLESS(number == 3,
                 "Number of vertices per face is expected to be 3");

    // Read in the identifier of the first face vertex.
    ERROR_UNLESS(lexer.get_integer(number),
                 "Error while trying to read in face vertex identifier");

    unsigned i1 = unsigned(number);

    ERROR_UNLESS(i1 < nv, "Face vertex identifier is invalid");

    // Read in the identifier of the second face vertex.
    ERROR_UNLESS(lexer.get_integer(number),
                 "Error while trying to read in face vertex identifier");

    unsigned i2 = unsigned(number);

    ERROR_UNLESS(i2 < nv, "Face vertex identifier is invalid");

    // Read in the identifier of the third face vertex.
    ERROR_UNLESS(lexer.get_integer(number),
                 "Error while trying to read in face vertex identifier");

    unsigned i3 = unsigned(number);

    ERROR_UNLESS(i3 < nv, "Face vertex identifier is invalid");

    // Store face data.
    const unsigned j = 3 * i;

    fset[j] = i1;
    fset[j + 1] = i2;
    fset[j + 2] = i3;
  }
}

} // namespace off
