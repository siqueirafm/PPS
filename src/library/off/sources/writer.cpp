/**
 * \file writer.cpp
 *
 * \brief Implementation of the  class Writer, which represents a file
 * writer for  writing out a  triangle surface mesh information  to an
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

#include "writer.hpp"

#include <exception.hpp>

#include <iomanip>
#include <iostream>

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

namespace fs = std::filesystem;

Writer::Writer(const fs::path &filepath) : _filepath(filepath) {
  const auto folder = _filepath.parent_path();

  if (fs::exists(folder)) {
    ERROR_UNLESS(fs::is_directory(folder),
                 "Path to the given file does not correspond to a folder");
  } else {
    ERROR_UNLESS(fs::create_directory(folder),
                 "Could not create folder " + folder.string());
  }

  const auto extension = filepath.extension();
  ERROR_UNLESS(extension == ".off", "Given file must have extension '.off'");
}

void Writer::write(unsigned nv, double *vset, unsigned nf, unsigned *fset) {
  _fs.open(_filepath.string().c_str(), std::ios::out | std::ios::binary);

  ERROR_UNLESS(_fs, "Could not open file " + _filepath.string());

  write_header(nv, nf);

  write_vertices(nv, vset);

  write_faces(nf, fset);

  _fs.close();
}

void Writer::write_header(unsigned nv, unsigned nf) {
  // Write first line
  _fs << "OFF" << std::endl;

  // Write the number of vertices, faces, and edges
  _fs << nv << " " << nf << " 0" << std::endl;
}

void Writer::write_vertices(unsigned nv, double *vset) {
  _fs << std::fixed << std::setprecision(18);

  for (unsigned i = 0; i < nv; i++) {

    const unsigned j = 3 * i;

    _fs << vset[j] << " " << vset[j + 1] << " " << vset[j + 2] << std::endl;
  }

  return;
}

void Writer::write_faces(unsigned nf, unsigned *fset) {
  for (unsigned i = 0; i < nf; i++) {

    const unsigned j = 3 * i;

    _fs << "3 " << fset[j] << " " << fset[j + 1] << " " << fset[j + 2]
        << std::endl;
  }

  return;
}

} // namespace off
