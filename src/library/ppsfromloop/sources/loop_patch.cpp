/**
 * \file loop_patch.cpp
 *
 * \brief  Implementation of  the class  LOOPatch, which  represents a
 * Loop subdivision surface patch.
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

#include "loop_patch.hpp" // LOOPatch

#include <exception.hpp> // ERROR_UNLESS

#include <cmath> // std::abs

namespace ppsfromloop {

LOOPatch::LOOPatch(EvalStruct *evs, unsigned evd, BCORDER bc, double *pp)
    : _evstruct(evs), _evdegree(evd), _bcorder(bc) {

  const unsigned int K = get_extraordinary_vertex_degree() + 6;

  _proj_pts = new double[K * 3];

  unsigned i;
  for (i = 0; i < K * 3; i++) {
    _proj_pts[i] = pp[i];
  }
}

LOOPatch::LOOPatch(const LOOPatch &patch) {
  if (&patch != this) {
    _evstruct = patch._evstruct;
    _evdegree = patch._evdegree;
    _bcorder = patch._bcorder;

    const unsigned int K = get_extraordinary_vertex_degree() + 6;

    _proj_pts = new double[K * 3];

    unsigned i;
    for (i = 0; i < K * 3; i++) {
      _proj_pts[i] = patch._proj_pts[i];
    }
  }
}

LOOPatch::LOOPatch(LOOPatch &&patch)
    : _evstruct(patch._evstruct), _evdegree(patch._evdegree),
      _bcorder(patch._bcorder), _proj_pts(std::move(patch._proj_pts)) {
  patch._proj_pts = nullptr;
  patch._evstruct = nullptr;
}

LOOPatch::~LOOPatch() {
  if (_proj_pts != 0) {
    delete[] _proj_pts;
  }
}

LOOPatch &LOOPatch::operator=(const LOOPatch &patch) {
  if (&patch != this) {
    _evstruct = patch._evstruct;
    _evdegree = patch._evdegree;
    _bcorder = patch._bcorder;

    const unsigned int K = get_extraordinary_vertex_degree() + 6;

    _proj_pts = new double[K * 3];

    unsigned i;
    for (i = 0; i < K * 3; i++) {
      _proj_pts[i] = patch._proj_pts[i];
    }
  }

  return *this;
}

EvalStruct *LOOPatch::get_surface_evaluator() const { return _evstruct; }

unsigned LOOPatch::get_extraordinary_vertex_degree() const { return _evdegree; }

LOOPatch::BCORDER LOOPatch::get_barycentric_coordinates_order() const {
  return _bcorder;
}

double *LOOPatch::get_projected_points() const { return _proj_pts; }

void LOOPatch::get_ith_projected_point(unsigned i, double &x, double &y,
                                       double &z) const {
  const unsigned int K = get_extraordinary_vertex_degree() + 6;

  ERROR_UNLESS(i < K, "Violated upper bound on degree of extraordinary vertex");

  x = _proj_pts[index(0, i, K)];
  y = _proj_pts[index(1, i, K)];
  z = _proj_pts[index(2, i, K)];
}

void LOOPatch::point(double u, double v, double &x, double &y,
                     double &z) const {
  if (std::abs(u) <= 1e-15) {
    u = 0;
  } else if (std::abs(1 - u) <= 1e-15) {
    u = 1;
  }

  if (std::abs(v) <= 1e-15) {
    v = 0;
  } else if (std::abs(1 - v) <= 1e-15) {
    v = 1;
  }

  double w = 1 - (u + v);

  if (std::abs(w) <= 1e-15) {
    w = 0;
  } else if (std::abs(1 - w) <= 1e-15) {
    w = 1;
  }

  double vv, ww;

  if (get_barycentric_coordinates_order() == UVW) {
    vv = v;
    ww = w;
  } else if (get_barycentric_coordinates_order() == VWU) {
    vv = w;
    ww = u;
  } else {
    vv = u;
    ww = v;
  }

  double s[3][3];
  get_surface_evaluator()->eval_surf(
      EvalStruct::EVAL_VALUE, get_projected_points(),
      get_extraordinary_vertex_degree(), vv, ww, s);

  x = s[0][0];
  y = s[0][1];
  z = s[0][2];
}

unsigned LOOPatch::index(unsigned i, unsigned j, unsigned n) {
  return j + (i * n);
}

} // namespace ppsfromloop
