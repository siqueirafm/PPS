#pragma once

/**
 * \file vertex_attribute.hpp
 *
 * \brief Definition  of class \c VertexAttribute,  which represents a
 * set of vertex attributes of the underlying triangle surface mesh of
 * a Parametric  Pseudo-Surface (PPS)  constructed from a  PN triangle
 * surface.
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

#include <bezier.hpp> // pps::Bezier

/**
 * \defgroup PPSFROMLOOPNameSpace Namespace ppsfromloop.
 * @{
 */

/**
 * \namespace ppsfromloop
 *
 * \brief  The  namespace  ppsfromloop  contains  the  definition  and
 * implementation   of   all   classes   representing   a   Parametric
 * Pseudo-Surface (PPS) that approximates the geometry of a given Loop
 * subdivision surface.
 */

namespace ppsfromloop {

using pps::Bezier;

/**
 * \class VertexAttribute
 *
 * \brief This class  represents a set of attributes  of a vertex of
 * the underlying triangle surface mesh  of a PPS constructed from a
 * Loop subdivision surface.
 */
class VertexAttribute {
public:
  /**
   * \fn VertexAttribute()
   *
   * \brief Creates an instance of this class.
   */
  VertexAttribute() = default;

  /**
   * \fn VertexAttribute( )
   *
   * \brief Creates an instance of this class.
   *
   * \param  patch  The  address   of  the  PPS  shape  function  (a
   * rectangular  B&eacute;zier patch )  associated with  the vertex
   * that owns this attribute.
   */
  explicit VertexAttribute(Bezier *patch);

  /**
   * \fn VertexAttribute( const VertexAttribute& a )
   *
   * \brief Creates an instance of this class from another instance.
   *
   * \param a A given instance of this class.
   */
  VertexAttribute(const VertexAttribute &a);

  /**
   * \fn VertexAttribute(VertexAttribute&& a )
   *
   * \brief Creates an instance of this class from another instance.
   *
   * \param a A given instance of this class.
   */
  VertexAttribute(VertexAttribute &&a);

  /**
   * \fn ~VertexAttribute()
   *
   * \brief Destroys an instance of this class.
   *
   */
  ~VertexAttribute();

  /**
   * \fn VertexAttribute& operator=( const VertexAttribute& a )
   *
   * \brief Overloads assignment operator.
   *
   * \param a A given instance of this class.
   */
  VertexAttribute &operator=(const VertexAttribute &a);

  /**
   * \fn Bezier* get_patch() const
   *
   * \brief  Returns  a  pointer  to  the PPS  shape  function  (  a
   * rectangular  B&eacute;zier patch)  associated  with the  vertex
   * that owns this attribute set.
   *
   * \return  A  pointer  to  the  rectangular  B&eacute;zier  patch
   * associated with the vertex that owns this attribute set.
   */
  Bezier *get_patch() const;

  /**
   * \fn void set_patch( Bezier* patch )
   *
   * \brief  Assigns  an  address  to  the  patch  pointer  of  this
   * attribute.
   *
   * \param  patch  The  address   of  the  PPS  shape  function  (a
   * rectangular  B&eacute;zier patch )  associated with  the vertex
   * that owns this attribute.
   */
  void set_patch(Bezier *patch);

private:
  /** Pointer to the shape function associated with the vertex that owns this
   * attribute. */
  Bezier *_patch = nullptr;
};

} // namespace ppsfromloop

/** @} */ // end of group class.
