#pragma once

/**
 * \file loop_patch.hpp
 *
 * \brief Definition  of class  \c LOOPatch,  which represents  a Loop
 * subdivision surface patch.
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

#include "loop_eval.hpp" // EvalStruct

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

/**
 * \class LOOPatch
 *
 * \brief This class represents a Loop subdivision surface patch.
 *
 */
class LOOPatch {
public:
  /**
   * \enum BCORDER
   * \brief  Constants   defining  the  order  of   the  barycentric
   * coordinates with  respect to  the affine frame  associated with
   * the base mesh face associated with this patch.
   */
  enum BCORDER { UVW = 0, VWU = 1, WUV = 2 };

  /**
   * \fn LOOPatch( EvalStruct* evs , unsigned evd , BCORDER bc , double* ppts )
   *
   * \brief Creates an instance of this class.
   *
   * \param  evs Pointer  to a  Loop subdivision  surface evaluator.
   * \param evd Degree of the  only extraordinary vertex (if any) of
   * the base  mesh face associated with  this face. If  there is no
   * extraordinary vertex, the default parameter value is 6.  \param
   * bc  Order of  the barycentric  coordinates associated  with the
   * affine  frame corresponding  to the  base mesh  face associated
   * with this patch.   \param pp Pointer to the  array of projected
   * points of this patch.
   */
  LOOPatch(EvalStruct *evs, unsigned evd, BCORDER bc, double *pp);

  /**
   * \fn LOOPatch( const LOOPatch& patch )
   *
   * \brief Creates an instance of this class from another instance.
   *
   * \param patch An instance of this class.
   */
  LOOPatch(const LOOPatch &patch);

  /**
   * \fn LOOPatch( LOOPatch&& patch )
   *
   * \brief Creates an instance of this class from another instance.
   *
   * \param patch An instance of this class.
   */
  LOOPatch(LOOPatch &&patch);

  /**
   * \fn ~LOOPatch()
   *
   * \brief Destroys an instance of this class.
   *
   */
  virtual ~LOOPatch();

  /**
   * \fn operator=( const LOOPatch& patch )
   *
   * \brief Overloads assignment operator.
   *
   * \param patch An instance of this class.
   *
   * \return A reference to the callee instance.
   */
  LOOPatch &operator=(const LOOPatch &patch);

  /**
   * \fn EvalStruct* get_surface_evaluator() const
   *
   * \brief  Returns a  pointer  to the  surface  evaluator of  this
   * patch.
   *
   * \return The pointer to the surface evaluator of this patch.
   */
  EvalStruct *get_surface_evaluator() const;

  /**
   * \fn unsigned get_extraordinary_vertex_degree() const
   *
   * \brief Returns  the degree of  the extraordinary vertex  of the
   * face that owns this attribute. If the face has no extraordinary
   * vertex, it returns the degree  of the first vertex of the face,
   * which should always be equal to 6.
   *
   * \return The degree of the extraordinary vertex of the base mesh
   * face associated with this patch.
   */
  unsigned get_extraordinary_vertex_degree() const;

  /**
   * \fn BCORDER get_barycentric_coordinates_order() const
   *
   * \brief  Returns  the   order  of  the  barycentric  coordinates
   * associated with the affine frame corresponding to the base mesh
   * face associated with patch.
   *
   * \return  The order  of the  barycentric  coordinates associated
   * with  the affine  frame  corresponding to  the  base mesh  face
   * associated with this patch.
   */
  BCORDER get_barycentric_coordinates_order() const;

  /**
   * \fn double* get_projected_points() const
   *
   * \brief Returns the projected points of this patch.
   *
   * \return A pointer to the projected points of this patch.
   */
  double *get_projected_points() const;

  /**
   * \fn void get_ith_projected_point( unsigned i , double& x , double& y ,
   * double& z ) const
   *
   * \brief Obtains the i-th projected  point of this patch.
   *
   * \param i Index of the point to be returned.
   * \param x First Cartesian coordinate of the i-th point.
   * \param y Second Cartesian coordinate of the i-th point.
   * \param z Third Cartesian coordinate of the i-th point.
   */
  void get_ith_projected_point(unsigned i, double &x, double &y,
                               double &z) const;

  /**
   * \fn virtual void point( double u , double v , double& x , double& y ,
   * double& z ) const
   *
   * \brief Computes a point on this patch.
   *
   * \param u First barycentric coordinate of a given parameter point.
   * \param v Second barycentric coordinate of a given parameter point.
   * \param x First Cartesian coordinate of the point on this patch.
   * \param y Second Cartesian coordinate of the point on this patch.
   * \param z Third Cartesian coordinate of the point on this patch.
   *
   */
  virtual void point(double u, double v, double &x, double &y, double &z) const;

  /**
   * \fn static unsigned index( unsigned i , unsigned j , unsigned n )
   *
   * \brief  Linearizes the  index (i,j)  of a  bidimensional matrix
   * with "n" columns.
   *
   * \param i Row number of the bidimensional matrix.
   * \param j Column number of the bidimensional matrix.
   * \param n Number of columns of the matrix.
   *
   * \return  Linear index corresponding  to the  given (row,column)
   * pair index.
   */
  static unsigned index(unsigned i, unsigned j, unsigned n);

private:
  // ---------------------------------------------------------------
  //
  // Private data members.
  //
  // ---------------------------------------------------------------

  /**
   * Pointer to the Loop subdivision surface evaluator.
   */
  EvalStruct *_evstruct = nullptr;

  /**
   * Degree of  the only extraordinary  vertex (if any) of  the base
   * mesh  face associated  with  this  patch. If  the  face has  no
   * extraordinary vertex,  then the value of the  attribute must be
   * 6.
   */
  unsigned _evdegree = 0;

  /**
   * Order of the barycentric coordinates associated with the affine
   * frame corresponding to the  base mesh face associated with this
   * patch.
   */
  BCORDER _bcorder = BCORDER::UVW;

  /**
   * A pointer to the array of projected points of this patch.
   */
  double *_proj_pts = nullptr;
};

} // namespace ppsfromloop

/** @} */ // end of group class.
