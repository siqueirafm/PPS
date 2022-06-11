#pragma once

/**
 * \file reader.hpp
 *
 * \brief  Definition of the  class Reader,  which represents  a file
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
 */

#include "lexer.hpp"

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
 * \class Reader
 *
 * \brief  This class  represents a  file  reader for  reading in  a
 * triangle surface mesh information from an OFF file.
 */
class Reader {
public:
  /**
   * \fn Reader( const std::filesystem::path& filepath )
   *
   * \brief Creates an instance of this class.
   *
   * \param filepath Full path to an OFF file describing a triangle mesh.
   */
  Reader(const std::filesystem::path &filepath);

  /**
   * \fn read( unsigned& nv , double*& vset , unsigned& nf , unsigned*& fset )
   *
   * \brief Reads  in the  topological and geometric  information of
   * the surface mesh described in the input file.
   *
   * \param nv A reference to the number of vertices of the surface.
   * \param  vset   A  reference  to   an  array  with   the  vertex
   * coordinates.
   * \param nf A reference to the number of faces of the surface.
   * \param  fset A  reference  to  an array  with  the face  vertex
   * identifiers.
   */
  void read(unsigned &nv, double *&vset, unsigned &nf, unsigned *&fset);

private:
  /**
   * \fn read_header(Lexer& lexer, unsigned& nv , unsigned& nf )
   *
   * \brief Reads in the number of vertices and faces of the surface
   * described by the input file.
   *
   * \param lexer A lexical analyzer for the input stream.
   * \param nv A reference to the number of vertices of the surface.
   * \param nf A reference to the number of faces of the surface.
   */
  void read_header(Lexer &lexer, unsigned &nv, unsigned &nf);

  /**
   * \fn read_vertices( Lexer& lexer, unsigned nv , double*& vset )
   *
   * \brief Reads in the vertex coordinates of the surface described
   * by the input file.
   *
   * \param lexer A lexical analyzer for the input stream.
   * \param nv The number of vertices of the surface.
   * \param vset A reference to an array of vertex coordinates.
   */
  void read_vertices(Lexer &lexer, unsigned nv, double *&vset);

  /**
   * \fn read_faces( Lexer& lexer, unsigned nf , unsigned nv , unsigned*& fset )
   *
   * \brief  Reads in  the face  vertex identifiers  of  the surface
   * described by the input file.
   *
   * \param lexer A lexical analyzer for the input stream.
   * \param nf The number of faces of the surface.
   * \param nv The number of vertices of the surface.
   * \param fset A reference to an array of face vertex identifiers.
   */
  void read_faces(Lexer &lexer, unsigned nf, unsigned nv, unsigned *&fset);

private:
  std::filesystem::path _filepath; ///< The input file name.
};

} // namespace off

/** @} */ // end of group class.
