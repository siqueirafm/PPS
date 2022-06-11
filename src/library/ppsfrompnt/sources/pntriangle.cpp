/**
 * \file pntriangle.cpp
 *
 * \brief Implementation  of the class PNTriangle,  which represents a
 * PN triangle.
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

#include "pntriangle.hpp"

#include <exception.hpp> // ERROR_UNLESS

#include <cmath> // std::abs, std::sqrt

namespace ppsfrompnt {

PNTriangle::PNTriangle(double p0[3], double p1[3], double p2[3], double n0[3],
                       double n1[3], double n2[3]) {
  // ---------------------------------------------------------------
  //
  // Compute control points b300, b030, and b003.
  //
  // ---------------------------------------------------------------

  /**
   * b300 = p0, b030 = p1, and b003 = p2
   */
  for (unsigned i = 0; i < 3; i++) {
    _b[0][i] = p0[i];
    _b[1][i] = p1[i];
    _b[2][i] = p2[i];
  }

  // ---------------------------------------------------------------
  //
  // Compute control points b210, b120, b201, b021, b102, and b012.
  //
  // ---------------------------------------------------------------

  /**
   * Compute the weight associated with the above control points.
   */

  /**
   * w210 = ( p1 - p0 ).inner( n0 ) ;
   */
  double w210 = 0;
  for (unsigned i = 0; i < 3; i++) {
    w210 += (p1[i] - p0[i]) * n0[i];
  }

  /**
   * w120 = ( p0 - p1 ).inner( n1 ) ;
   */
  double w120 = 0;
  for (unsigned i = 0; i < 3; i++) {
    w120 += (p0[i] - p1[i]) * n1[i];
  }

  /**
   * w201 = ( p2 - p0 ).inner( n0 ) ;
   */
  double w201 = 0;
  for (unsigned i = 0; i < 3; i++) {
    w201 += (p2[i] - p0[i]) * n0[i];
  }

  /**
   * w021 = ( p2 - p1 ).inner( n1 ) ;
   */
  double w021 = 0;
  for (unsigned i = 0; i < 3; i++) {
    w021 += (p2[i] - p1[i]) * n1[i];
  }

  /**
   * w102 = ( p0 - p2 ).inner( n2 ) ;
   */
  double w102 = 0;
  for (unsigned i = 0; i < 3; i++) {
    w102 += (p0[i] - p2[i]) * n2[i];
  }

  /**
   *  w012 = ( p1 - p2 ).inner( n2 ) ;
   */
  double w012 = 0;
  for (unsigned i = 0; i < 3; i++) {
    w012 += (p1[i] - p2[i]) * n2[i];
  }

  /**
   * b210 = ( 1. / 3. ) * ( ( 2 * p0 ) + p1 - ( w210 * n0 ) ) ;
   */
  for (unsigned i = 0; i < 3; i++) {
    _b[3][i] = (2 * p0[i]) + p1[i] - (w210 * n0[i]);
    _b[3][i] /= 3.;
  }

  /**
   * b120 = ( 1. / 3. ) * ( ( 2 * p1 ) + p0 - ( w120 * n1 ) ) ;
   */
  for (unsigned i = 0; i < 3; i++) {
    _b[4][i] = (2 * p1[i]) + p0[i] - (w120 * n1[i]);
    _b[4][i] /= 3.;
  }

  /**
   * b201 = ( 1. / 3. ) * ( ( 2 * p0 ) + p2 - ( w201 * n0 ) ) ;
   */
  for (unsigned i = 0; i < 3; i++) {
    _b[5][i] = (2 * p0[i]) + p2[i] - (w201 * n0[i]);
    _b[5][i] /= 3.;
  }

  /**
   * b021 = ( 1. / 3. ) * ( ( 2 * p1 ) + p2 - ( w021 * n1 ) ) ;
   */
  for (unsigned i = 0; i < 3; i++) {
    _b[6][i] = (2 * p1[i]) + p2[i] - (w021 * n1[i]);
    _b[6][i] /= 3.;
  }

  /**
   * b102 = ( 1. / 3. ) * ( ( 2 * p2 ) + p0 - ( w102 * n2 ) ) ;
   */
  for (unsigned i = 0; i < 3; i++) {
    _b[7][i] = (2 * p2[i]) + p0[i] - (w102 * n2[i]);
    _b[7][i] /= 3.;
  }

  /**
   * b012 = ( 1. / 3. ) * ( ( 2 * p2 ) + p1 - ( w012 * n2 ) ) ;
   */
  for (unsigned i = 0; i < 3; i++) {
    _b[8][i] = (2 * p2[i]) + p1[i] - (w012 * n2[i]);
    _b[8][i] /= 3.;
  }

  // ---------------------------------------------------------------
  //
  // Compute control point b111.
  //
  // ---------------------------------------------------------------

  /**
   * b111 = e + ( ( 0.5 * ( e - v ) ) )
   *
   * where
   *
   * e = ( b210 + b120 + b021 + b012 + b102 + b201 ) / 6
   *
   * and
   *
   * v = ( p0 + p1 + p2 ) / 3
   *
   */
  double e[3];
  double v[3];
  for (unsigned i = 0; i < 3; i++) {
    e[i] = 0.0;
    for (unsigned j = 3; j < 9; j++) {
      e[i] += _b[j][i];
    }

    e[i] /= 6.;

    v[i] = (p0[i] + p1[i] + p2[i]) / 3.;

    _b[9][i] = e[i] + (0.5 * (e[i] - v[i]));
  }

  // ---------------------------------------------------------------
  //
  // Compute unit normals n200, n020, and n002.
  //
  // ---------------------------------------------------------------

  /**
   * n200 = n0, n020 = n1, and n002 = n2.
   */
  for (unsigned i = 0; i < 3; i++) {
    _n[0][i] = n0[i];
    _n[1][i] = n1[i];
    _n[2][i] = n2[i];
  }

  // ---------------------------------------------------------------
  //
  // Compute unit normals n110, n101, and n011.
  //
  // ---------------------------------------------------------------

  /**
   * Compute the weights associated with each unit normal.
   */

  /**
   * n110 = mv1 - ( ( 2 * w110 ) / d110 ) * ( mv0 ) ;
   */

  /**
   * mv0 = p1 - p0
   * mv1 = n0 + n1
   */
  double mv0[3];
  double mv1[3];
  for (unsigned i = 0; i < 3; i++) {
    mv0[i] = p1[i] - p0[i];
    mv1[i] = n0[i] + n1[i];
  }

  /**
   * d110 = mv0.inner( mv0 )
   * w110 = mv0.inner( mv1 )
   */
  double d110 = 0;
  double w110 = 0;
  for (unsigned i = 0; i < 3; i++) {
    d110 += mv0[i] * mv0[i];
    w110 += mv0[i] * mv1[i];
  }

  for (unsigned i = 0; i < 3; i++) {
    _n[3][i] = mv1[i] - ((2 * w110) / d110) * (mv0[i]);
  }

  /**
   * n101 = mv1 - ( ( 2 * w101 ) / d101 ) * ( mv0 ) ;
   */

  /**
   * mv0 = p2 - p0
   * mv1 = n0 + n2
   */
  for (unsigned i = 0; i < 3; i++) {
    mv0[i] = p2[i] - p0[i];
    mv1[i] = n0[i] + n2[i];
  }

  /**
   * d101 = mv0.inner( mv0 )
   * w101 = mv0.inner( mv1 )
   */
  double d101 = 0;
  double w101 = 0;
  for (unsigned i = 0; i < 3; i++) {
    d101 += mv0[i] * mv0[i];
    w101 += mv0[i] * mv1[i];
  }

  for (unsigned i = 0; i < 3; i++) {
    _n[4][i] = mv1[i] - ((2 * w101) / d101) * (mv0[i]);
  }

  /**
   * n011 = mv1 - ( ( 2 * w011 ) / d011 ) * ( mv0 ) ;
   */

  /**
   * mv0 = p2 - p1
   * mv1 = n1 + n2
   */
  for (unsigned i = 0; i < 3; i++) {
    mv0[i] = p2[i] - p1[i];
    mv1[i] = n1[i] + n2[i];
  }

  /**
   * d011 = mv0.inner( mv0 )
   * w011 = mv0.inner( mv1 )
   */
  double d011 = 0;
  double w011 = 0;
  for (unsigned i = 0; i < 3; i++) {
    d011 += mv0[i] * mv0[i];
    w011 += mv0[i] * mv1[i];
  }

  for (unsigned i = 0; i < 3; i++) {
    _n[5][i] = mv1[i] - ((2 * w011) / d011) * (mv0[i]);
  }

  /**
   * Normalize the vectors n110, n101, and n011.
   */
  double l[3];
  for (unsigned i = 3; i < 6; i++) {
    l[i - 3] = 0.0;

    for (unsigned j = 0; j < 3; j++) {
      l[i - 3] += _n[i][j] * _n[i][j];
    }

    double sqrtl = std::sqrt(l[i - 3]);

    for (unsigned j = 0; j < 3; j++) {
      _n[i][j] /= sqrtl;
    }
  }

  return;
}

void PNTriangle::point(double u, double v, double &x, double &y,
                       double &z) const {
  if (fabs(u) <= 1e-15) {
    u = 0;
  } else if (fabs(1 - u) <= 1e-15) {
    u = 1;
  }

  if (fabs(v) <= 1e-15) {
    v = 0;
  } else if (fabs(1 - v) <= 1e-15) {
    v = 1;
  }

  double w = 1 - (u + v);

  if (fabs(w) <= 1e-15) {
    w = 0;
  } else if (fabs(1 - w) <= 1e-15) {
    w = 1;
  }

  double u2 = u * u;
  double v2 = v * v;
  double w2 = w * w;

  double u3 = u * u2;
  double v3 = v * v2;
  double w3 = w * w2;

  double pt[3];
  for (unsigned i = 0; i < 3; i++) {
    pt[i] = u3 * _b[0][i] + v3 * _b[1][i] + w3 * _b[2][i] +
            (3 * u2 * v) * _b[3][i] + (3 * u * v2) * _b[4][i] +
            (3 * u2 * w) * _b[5][i] + (3 * v2 * w) * _b[6][i] +
            (3 * u * w2) * _b[7][i] + (3 * v * w2) * _b[8][i] +
            (6 * u * v * w) * _b[9][i];
  }

  x = pt[0];
  y = pt[1];
  z = pt[2];

  return;
}

void PNTriangle::normal(double u, double v, double &x, double &y,
                        double &z) const {

  if (fabs(u) <= 1e-15) {
    u = 0;
  } else if (fabs(1 - u) <= 1e-15) {
    u = 1;
  }

  if (fabs(v) <= 1e-15) {
    v = 0;
  } else if (fabs(1 - v) <= 1e-15) {
    v = 1;
  }

  double w = 1 - (u + v);

  if (fabs(w) <= 1e-15) {
    w = 0;
  } else if (fabs(1 - w) <= 1e-15) {
    w = 1;
  }

  ERROR_UNLESS((u >= 0) && (u <= 1) && (v >= 0) && (v <= 1) && (w >= 0) &&
                   (w <= 1),
               "Barycentric coordinates are invalid");

  double u2 = u * u;
  double v2 = v * v;
  double w2 = w * w;

  double pt[3];
  for (unsigned i = 0; i < 3; i++) {
    pt[i] = u2 * _n[0][i] + v2 * _n[1][i] + w2 * _n[2][i] + (u * v) * _n[3][i] +
            (u * w) * _n[4][i] + (v * w) * _n[5][i];
  }

  x = pt[0];
  y = pt[1];
  z = pt[2];

  double l = std::sqrt((x * x) + (y * y) + (z * z));

  ERROR_UNLESS(l > 0, "Normal vector cannot be the null vector");

  x /= l;
  y /= l;
  z /= l;

  return;
}

} // namespace ppsfrompnt
