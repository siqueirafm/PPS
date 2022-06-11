#pragma once

/**
 * \file face_attribute.hpp
 *
 * \brief Definition of class \c FaceAttribute, which represents a set
 * of face  attributes of  the underlying triangle  surface mesh  of a
 * Parametric Pseudo-Surface (PPS) constructed from a Loop subdivision
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

#include "loop_patch.hpp" // LOOPatch

/**
 * \defgroup PPSFROMLOOPNameSpace Namespace ppsfromloop.
 * @{
 */

/**
 * \namespace ppsfromloop
 *
 * \brief  The  namespace   ppsfromloop  contains  the  definition  and
 * implementation   of   all    classes   representing   a   Parametric
 * Pseudo-Surface (PPS) that approximates  the geometry of a given Loop
 * subdivision surface.
 */

namespace ppsfromloop {

/**
 * \class FaceAttribute
 *
 * \brief This class represents a set of attributes of a face of the
 * underlying triangle surface mesh of a PPS constructed from a Loop
 * subdivision surface.
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
   * \fn FaceAttribute( LOOPatch* patch )
   *
   * \brief Creates an instance of this class.
   *
   * \param patch  A pointer to  the Loop subdivision  surface patch
   * associated with the face that owns this attribute.
   */
  FaceAttribute(LOOPatch *patch);

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
   * \fn LOOPatch* get_patch() const
   *
   * \brief  Returns  a  pointer  to the  triangular  surface  patch
   * associated with this face.
   *
   * \return A  pointer to  the triangular surface  patch associated
   * with this face.
   */
  LOOPatch *get_patch() const;

  /**
   * \fn void set_patch( LOOPatch* patch )
   *
   * \brief  Assigns an  address to  the pointer  to  the triangular
   * surface patch associated with this face.
   *
   * \param patch The address  of the Loop subdivision surface patch
   * associated with this face.
   */
  void set_patch(LOOPatch *patch);

private:
  /** Pointer to the Loop subdivision surface patch associated with this face.
   */
  LOOPatch *_patch;
};

} // namespace ppsfromloop

/** @} */ // end of group class.
