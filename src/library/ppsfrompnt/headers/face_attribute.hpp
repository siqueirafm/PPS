#pragma once

/**
 * \file face_attribute.hpp
 *
 * \brief Definition  of the  class FaceAttribute, which  represents a
 * set of face attributes of the underlying triangle surface mesh of a
 * Parametric  Pseudo-Surface  (PPS) constructed  from  a PN  triangle
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
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#include "pntriangle.hpp" // PNTriangle

/**
 * \defgroup PPSFROMPNTNameSpace Namespace ppsfrompnt.
 * @{
 */

/**
 * \namespace ppsfrompnt
 *
 * \brief  The  namespace  ppsfrompnt   contains  the  definition  and
 * implementation   of   all   classes   representing   a   Parametric
 * Pseudo-Surface (PPS) that  approximates the geometry of  a given PN
 * triangle surface.
 */

namespace ppsfrompnt {

/**
 * \class FaceAttribute
 *
 * \brief This class represents a set of attributes of a face of the
 * underlying triangle surface  mesh of a PPS constructed  from a PN
 * triangle surface.
 */
class FaceAttribute {
public:
  /**
   * \fn FaceAttribute()
   *
   * \brief Creates an instance of this class.
   */
  FaceAttribute();

  /**
   * \fn FaceAttribute( PNTriangle* patch )
   *
   * \brief Creates an instance of this class.
   *
   * \param patch A  pointer to the PN triangle  associated with the
   * face that owns this attribute.
   */
  FaceAttribute(PNTriangle *patch);

  /**
   * \fn FaceAttribute( const FaceAttribute& a )
   *
   * \brief Creates an instance of this class from another instance.
   *
   * \param a A given instance of this class.
   */
  FaceAttribute(const FaceAttribute &a);

  /**
   * \fn FaceAttribute( FaceAttribute&& a )
   *
   * \brief Creates an instance of this class from another instance.
   *
   * \param a A given instance of this class.
   */
  FaceAttribute(FaceAttribute &&a);

  /**
   * \fn ~FaceAttribute()
   *
   * \brief Destroys an instance of this class.
   *
   */
  ~FaceAttribute();

  /**
   * \fn FaceAttribute& operator=( const FaceAttribute& a )
   *
   * \brief Overloads assignment operator.
   *
   * \param a A given instance of this class.
   */
  FaceAttribute &operator=(const FaceAttribute &a);

  /**
   * \fn PNTriangle* get_patch() const
   *
   * \brief  Returns  a  pointer  to the  triangular  surface  patch
   * associated with this face.
   *
   * \return A  pointer to  the triangular surface  patch associated
   * with this face.
   */
  PNTriangle *get_patch() const;

  /**
   * \fn void set_patch( PNTriangle* patch )
   *
   * \brief  Assigns an  address to  the pointer  to  the triangular
   * surface patch associated with this face.
   *
   * \param  patch The address  of the  PN triangle  associated with
   * this face.
   */
  void set_patch(PNTriangle *patch);

private:
  /** Pointer to the PN triangle associated with this face. */
  PNTriangle *_patch;
};

} // namespace ppsfrompnt

/** @} */ // end of group class.
