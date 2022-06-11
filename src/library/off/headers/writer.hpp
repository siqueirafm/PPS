#pragma once

/**
 * \file writer.hpp
 *
 * \brief  Definition of  the class  Writer, which  represents  a file
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

#include <filesystem>
#include <fstream>
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
 * \class Writer
 *
 * \brief  This class  represents a  file writer  for writing  out a
 * triangle surface mesh information to an OFF file.
 */
class Writer {
public:
  /**
   * \fn Writer( const std::filesystem::path& filepath )
   *
   * \brief Creates an instance of this class.
   *
   * \param fn The name of an OFF input file.
   */
  Writer(const std::filesystem::path &filepath);

  /**
   * \fn void write( unsigned nv , double* vset , unsigned nf , unsigned* fset )
   *
   * \brief Writes the vertex and face information of a surface mesh
   * to an OFF file.
   *
   * \param nv The number of vertices of the mesh.
   * \param vset The Cartesian coordinates of the mesh vertices.
   * \param nf The number of faces of the mesh.
   * \param fset The set of vertex identifiers of the mesh faces.
   */
  void write(unsigned nv, double *vset, unsigned nf, unsigned *fset);

private:
  /**
   * \fn void write_header( unsigned nv , unsigned nf )
   *
   * \brief Writes the file header information.
   *
   * \param nv The number of vertices of the mesh.
   * \param nf The number of faces of the mesh.
   */
  void write_header(unsigned nv, unsigned nf);

  /**
   * \fn void write_vertices( unsigned nv , double* vset )
   *
   * \brief Writes the file header information.
   *
   * \param nv The number of vertices of the mesh.
   * \param vset The Cartesian coordinates of the mesh vertices.
   */
  void write_vertices(unsigned nv, double *vset);

  /**
   * \fn void write_faces( unsigned nf , unsigned* fset )
   *
   * \brief Writes the mesh face information.
   *
   * \param nf The number of faces of the mesh.
   * \param fset The set of vertex identifiers of the mesh faces.
   */
  void write_faces(unsigned nf, unsigned *fset);

  // ---------------------------------------------------------------
  //
  // Private data members
  //
  // ---------------------------------------------------------------

  std::filesystem::path _filepath; ///< Full path to the OFF file to which mesh
                                   ///< data will be written.
  std::fstream _fs;                ///< The output file stream.
};

} // namespace off

/** @} */ // end of group class.
