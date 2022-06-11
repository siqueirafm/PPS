#pragma once

/**
 * \file pntriangle.hpp
 *
 * \brief Definition  of the class  PNTriangle, which represents  a PN
 * triangle.
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

/**
 * \defgroup PPSFROMPNTNameSpace Namespace ppsfrompnt.
 * @{
 */

/**
 * \namespace ppsfrompnt
 *
 * \brief  The   namespace  ppsfrompnt  contains   the  definition  and
 * implementation   of   all    classes   representing   a   Parametric
 * Pseudo-Surface (PPS)  that approximates the  geometry of a  given PN
 * triangle surface.
 */

namespace ppsfrompnt {

/**
 * \class PNTriangle
 *
 * \brief This class represents a PN triangle.
 *
 */
class PNTriangle {
public:
  /**
   * \fn PNTriangle( double p0[3] , double p1[3] , double p2[3] , double n0[3] ,
   * double n1[3] , double n2[3] )
   *
   * \brief Creates  an instance of  this class.
   *
   * \param  p0 Cartesian  coordinates of  the first  vertex  of the
   * triangle defining the PN triangle.
   * \param  p1 Cartesian coordinates  of the  second vertex  of the
   * triangle defining the PN triangle.
   * \param  p2 Cartesian  coordinates of  the third  vertex  of the
   * triangle defining the PN triangle.
   * \param n0 A unit normal at the first vertex.
   * \param n1 A unit normal at the second vertex.
   * \param n2 A unit normal at the third vertex.
   */
  PNTriangle(double p0[3], double p1[3], double p2[3], double n0[3],
             double n1[3], double n2[3]);

  /**
   * \fn void point( double u , double v , double& x , double& y , double& z )
   * const
   *
   * \brief Computes a point on this patch.
   *
   * \param u First barycentric coordinate of a given parameter point.
   * \param v Second barycentric coordinate of a given parameter point.
   * \param x First Cartesian coordinate of the point on this PN triangle.
   * \param y Second Cartesian coordinate of the point on this PN triangle.
   * \param z Third Cartesian coordinate of the point on this PN triangle.
   *
   */
  void point(double u, double v, double &x, double &y, double &z) const;

  /**
   * \fn void normal( double u , double v , double& x , double& y , double& z )
   * const
   *
   * \brief Computes the unit normal  of this PN triangle at a given
   * point.
   *
   * \param u First barycentric coordinate of a given parameter point.
   * \param v Second barycentric coordinate of a given parameter point.
   * \param x First component of the unit normal.
   * \param y Second component of the unit normal.
   * \param z Third component of the unit normal.
   */
  void normal(double u, double v, double &x, double &y, double &z) const;

private:
  /**
   * The "b( i , j , k )" coefficients of this PN triangle.
   *
   * _b[ 0 ] --> b300
   * _b[ 1 ] --> b030
   * _b[ 2 ] --> b003
   * _b[ 3 ] --> b210
   * _b[ 4 ] --> b120
   * _b[ 5 ] --> b201
   * _b[ 6 ] --> b021
   * _b[ 7 ] --> b102
   * _b[ 8 ] --> b012
   * _b[ 9 ] --> b111
   *
   */
  double _b[10][3];

  /**
   * The "n( i , j , k )" coefficients for computing the unit normal
   * of this PN triangle at a given point. These coefficients define
   * a quadratic interpolant.
   *
   * _n[ 0 ] --> n200
   * _n[ 1 ] --> n020
   * _n[ 2 ] --> n002
   * _n[ 3 ] --> n110
   * _n[ 4 ] --> n101
   * _n[ 5 ] --> n011
   *
   */
  double _n[6][3];
};

} // namespace ppsfrompnt

/** @} */ // end of group class.
