#pragma once

/**
 * \file face.hpp
 *
 * \brief Definition of the class Face, which represents a face (i.e.,
 * a  triangle)  from  a  surface  mesh represented  by  a  DCEL  data
 * structure.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 2.0
 * \date January 2010
 */

#include "halfedge.hpp" // Halfedge

/**
 * \defgroup DCELNameSpace Namespace dcel.
 * @{
 */

/**
 * \namespace dcel
 *
 * \brief   The   namespace   dcel   contains   the   definition   and
 * implementation of all classes  of the <strong>Doubly Connected Edge
 * List (DCEL)</strong> data structure, which can be used to represent
 * surface meshes with empty boundary.
 *
 */

namespace dcel {

/**
 * \class Face
 *
 * \brief This  class represents  a face (i.e.,  a triangle)  from a
 * surface mesh represented by the DCEL data structure.
 */
template <typename VAttrib, typename FAttrib, typename EAttrib,
          typename HAttrib>
class Face {
public:
  // ---------------------------------------------------------------
  //
  // Type definitions
  //
  // ---------------------------------------------------------------

  /**
   * \using Halfedge
   *
   * \brief  Defines  Halfedge   as  an  alias  for  dcel::Halfedge<
   * VAttrib, FAttrib , EAttrib , HAttrib >.
   */
  using Halfedge = dcel::Halfedge<VAttrib, FAttrib, EAttrib, HAttrib>;

  // ---------------------------------------------------------------
  //
  // Public methods
  //
  // ---------------------------------------------------------------

  /**
   * \fn Face( Halfedge* h )
   *
   * \brief Creates an instance of this class.
   *
   * \param h A pointer to the first halfedge of this face.
   */
  Face(Halfedge *h) : _halfedge(h) {}

  /**
   * \fn ~Face()
   *
   * \brief Destroys an instance of this class.
   */
  ~Face() { set_halfedge(0); }

  /**
   * \fn Halfedge* get_halfedge() const
   *
   * \brief Returns a pointer to the first half-edge of this face.
   *
   * \return A pointer to the first half-edge of this face.
   */
  Halfedge *get_halfedge() const { return _halfedge; }

  /**
   * \fn void set_halfedge( Halfedge* h )
   *
   * \brief Assigns an address to the pointer to the first half-edge
   * of this face.
   *
   * \param h The address of the first half-edge of this face.
   */
  void set_halfedge(Halfedge *h) { _halfedge = h; }

  /**
   * \fn FAttrib& get_attributes()
   *
   * \brief Returns the set of attributes associated with this face.
   *
   * \return A reference to the set of attributes associated with this
   * face.
   */
  FAttrib &get_attributes() { return _attributes; }

private:
  // ---------------------------------------------------------------
  //
  // Private data members
  //
  // ---------------------------------------------------------------

  Halfedge *_halfedge = nullptr; ///< Pointer to the first half-edge of the
                                 ///< half-edge cycle of this face.
  FAttrib _attributes; ///< Set of attributes associated with this face.
};

} // namespace dcel

/** @} */ // end of group class.
