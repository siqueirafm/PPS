/**
 * \file bezier.cpp
 *
 * \brief  Implementation  of  the  class  Bezier  that  represents  a
 * rectangular B&eacute;zier patch in 3D.
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

#include "bezier.hpp" // Bezier

#include "ludcmp.hpp" // LUdcmp

#include <exception.hpp> // ERROR_UNLESS

#include <cmath> // std::abs

namespace pps {

Bezier::Bezier(double *param_pts, double *image_pts, unsigned np, unsigned m,
               unsigned n, double rx, double ry, double sx, double sy) {
  /**
   * Check consistency of the input parameters.
   */
  ERROR_UNLESS(m > 0, "Patch degree must be a positive number");
  ERROR_UNLESS(n > 0, "Patch degree must be a positive number");
  ERROR_UNLESS(rx < sx, "Invalid affine frame");
  ERROR_UNLESS(ry < sy, "Invalid affine frame");
  ERROR_UNLESS((2 * (n + 1) * (m + 1)) <= np,
               "Insufficient number of sample points");

  /**
   * Initialize private attributes.
   */

  set_bidegree_1st_index(m);
  set_bidegree_2nd_index(n);

  set_aff_1st_point_x_coord(rx);
  set_aff_1st_point_y_coord(ry);

  set_aff_2nd_point_x_coord(sx);
  set_aff_2nd_point_y_coord(sy);

  /**
   * Allocate memory for the control points.
   */
  unsigned ncp =
      (get_bidegree_1st_index() + 1) * (get_bidegree_2nd_index() + 1);

  _ctrl_pts = new double *[ncp];

  for (unsigned i = 0; i < ncp; i++) {
    _ctrl_pts[i] = new double[3];
  }

  /**
   * Compute a (np x ncp)  matrix "A" with all Bernstein polynomials
   * of order "_m  + _n" at the parameter points  given in the input
   * vector "param_pts".
   */

  double **A = new double *[np];

  for (unsigned i = 0; i < np; i++) {
    A[i] = new double[ncp];
  }

  comp_bpoly_matrix(A, param_pts, np);

  /**
   * Compute the  product "At x A",  where "At" is  the transpose of
   * "A".
   */
  double **AtA = new double *[ncp];

  for (unsigned i = 0; i < ncp; i++) {
    AtA[i] = new double[ncp];
  }

  comp_matrix_ata(A, np, ncp, AtA);

  /**
   * Compute the  product "At x B",  where "At" is  the transpose of
   * "A".
   */
  double **AtB = new double *[ncp];

  for (unsigned i = 0; i < ncp; i++) {
    AtB[i] = new double[3];
  }

  comp_matrix_atb(A, image_pts, np, ncp, AtB);

  /**
   * Solve the linear system "AtA Y = AtB[i]", where "AtB[i]" is the
   * i-th column  of "AtB",  for i=0,1,2. This  is done by  using LU
   * decomposition.
   */

  /**
   * Compute the LU decomposition of AtA.
   */
  LUdcmp LU(AtA, ncp);

  /**
   * Allocate  memory for the  solution vector,  "Y", of  the linear
   * system "AtA Y = AtB[i].
   *
   */

  double *Y = new double[ncp];

  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < ncp; j++) {
      Y[j] = AtB[j][i];
    }

    /**
     * Solve the linear system AtA Y = AtB[i].
     */
    LU.solve(Y, ncp);

    /**
     * Store  the  solution  in   the  i-th  column  of  the  matrix
     * "_ctrl_pts[j]".
     */
    for (unsigned j = 0; j < ncp; j++) {
      _ctrl_pts[j][i] = Y[j];
    }
  }

  /**
   * Release memory
   */
  if (A != 0) {
    for (unsigned i = 0; i < np; i++) {
      if (A[i] != 0) {
        delete[] A[i];
      }
    }

    delete A;
  }

  if (AtA != 0) {
    for (unsigned i = 0; i < ncp; i++) {
      if (AtA[i] != 0) {
        delete[] AtA[i];
      }
    }

    delete AtA;
  }

  if (AtB != 0) {
    for (unsigned i = 0; i < ncp; i++) {
      if (AtB[i] != 0) {
        delete[] AtB[i];
      }
    }

    delete AtB;
  }

  delete[] Y;

  return;
}

Bezier::Bezier(const Bezier &bz) {
  /**
   * Initialize private attributes.
   */

  set_bidegree_1st_index(bz.get_bidegree_1st_index());
  set_bidegree_2nd_index(bz.get_bidegree_2nd_index());

  set_aff_1st_point_x_coord(bz.get_aff_1st_point_x_coord());
  set_aff_1st_point_y_coord(bz.get_aff_1st_point_y_coord());

  set_aff_2nd_point_x_coord(bz.get_aff_2nd_point_x_coord());
  set_aff_2nd_point_y_coord(bz.get_aff_2nd_point_y_coord());

  /**
   * Allocate memory for the control points.
   */
  unsigned ncp =
      (get_bidegree_1st_index() + 1) * (get_bidegree_2nd_index() + 1);

  _ctrl_pts = new double
      *[(get_bidegree_1st_index() + 1) * (get_bidegree_2nd_index() + 1)];

  for (unsigned i = 0; i < ncp; i++) {
    _ctrl_pts[i] = new double[3];

    for (unsigned j = 0; j < 3; j++) {
      _ctrl_pts[i][j] = bz._ctrl_pts[i][j];
    }
  }

  return;
}

Bezier::~Bezier() {
  if (_ctrl_pts != 0) {

    unsigned ncp =
        (get_bidegree_1st_index() + 1) * (get_bidegree_2nd_index() + 1);

    for (unsigned i = 0; i < ncp; i++) {
      if (_ctrl_pts[i] != 0) {
        delete[] _ctrl_pts[i];
      }
    }

    delete[] _ctrl_pts;
  }

  return;
}

void Bezier::b(unsigned i, unsigned j, double &x, double &y, double &z) const {
  ERROR_UNLESS((i <= get_bidegree_1st_index()) &&
                   (j <= get_bidegree_2nd_index()),
               "Invalid index of control point");

  unsigned l = index(j, i);

  x = _ctrl_pts[l][0];
  y = _ctrl_pts[l][1];
  z = _ctrl_pts[l][2];

  return;
}

void Bezier::point(double u, double v, double &x, double &y, double &z) const {
  /**
   * Map the point to the affine frame [0,1].
   */
  double rx = get_aff_1st_point_x_coord();
  double ry = get_aff_1st_point_y_coord();
  double sx = get_aff_2nd_point_x_coord();
  double sy = get_aff_2nd_point_y_coord();

  /**
   * Compute all Bernstein polynomials of degree \var{Bezier::_m}.
   */

  double uu = (u - rx) / (sx - rx);
  double vv = (v - ry) / (sy - ry);

  if (std::abs(uu) <= 1e-15) {
    uu = 0;
  } else if (std::abs(1 - uu) <= 1e-15) {
    uu = 1;
  }

  if (std::abs(vv) <= 1e-15) {
    vv = 0;
  } else if (std::abs(1 - vv) <= 1e-15) {
    vv = 1;
  }

  unsigned m = get_bidegree_1st_index();
  unsigned n = get_bidegree_2nd_index();

  std::vector<double> bu(m + 1);
  all_bernstein(m, uu, bu);

  /**
   * Compute all Bernstein polynomials of degree \var{Bezier::_n}.
   */
  std::vector<double> bv(n + 1);
  all_bernstein(n, vv, bv);

  /**
   * Compute the image point, (x,y,z), of (u,v) on this patch.
   */
  x = 0;
  y = 0;
  z = 0;
  for (unsigned j = 0; j <= n; j++) {
    for (unsigned i = 0; i <= m; i++) {
      double buv = bu[i] * bv[j];
      double xaux;
      double yaux;
      double zaux;
      b(i, j, xaux, yaux, zaux);

      x += buv * xaux;
      y += buv * yaux;
      z += buv * zaux;
    }
  }

  return;
}

double Bezier::bernstein(unsigned n, unsigned i, double u) const {
  ERROR_UNLESS(i <= n, "Invalid index of Bernstein polynomial");

  std::vector<double> temp(n + 1);

  for (unsigned j = 0; j <= n; j++) {
    temp[j] = 0;
  }

  temp[n - i] = 1;

  double u1 = 1 - u;

  for (unsigned k = 1; k <= n; k++) {
    for (unsigned j = n; j >= k; j--) {
      temp[j] = (u1 * temp[j]) + (u * temp[j - 1]);
    }
  }

  return temp[n];
}

void Bezier::all_bernstein(unsigned n, double u, std::vector<double> &b) const {
  b[0] = 1;
  double u1 = 1 - u;

  for (unsigned j = 1; j <= n; j++) {
    double saved = 0;
    for (unsigned k = 0; k < j; k++) {
      double temp = b[k];
      b[k] = saved + (u1 * temp);
      saved = u * temp;
    }

    b[j] = saved;
  }

  return;
}

void Bezier::comp_bpoly_matrix(double **&a, double *param_pts,
                               unsigned np) const {
  /**
   * Compute the value of the Bernstein polynomials at the parameter
   * points.
   */
  double rx = get_aff_1st_point_x_coord();
  double ry = get_aff_1st_point_y_coord();
  double sx = get_aff_2nd_point_x_coord();
  double sy = get_aff_2nd_point_y_coord();

  for (unsigned i = 0; i < np; i++) {
    double x = (param_pts[(i << 1)] - rx) / (sx - rx);

    double y = (param_pts[(i << 1) + 1] - ry) / (sy - ry);

    unsigned m = get_bidegree_1st_index();
    unsigned n = get_bidegree_2nd_index();

    std::vector<double> bu(m + 1);
    std::vector<double> bv(n + 1);

    all_bernstein(m, x, bu);
    all_bernstein(n, y, bv);

    for (unsigned k = 0; k <= n; ++k) {
      for (unsigned j = 0; j <= m; ++j) {
        unsigned l = index(k, j);
        a[i][l] = bu[j] * bv[k];
      }
    }
  }

  return;
}

void Bezier::comp_matrix_ata(double **a, unsigned n, unsigned p,
                             double **&ata) const {
  for (unsigned i = 0; i < p; i++) {
    for (unsigned k = 0; k < p; k++) {
      ata[i][k] = 0;
      for (unsigned j = 0; j < n; j++) {
        ata[i][k] += a[j][i] * a[j][k];
      }
    }
  }

  return;
}

void Bezier::comp_matrix_atb(double **a, double *b, unsigned n, unsigned p,
                             double **&atb) const {
  /**
   * Compute the  product \f$a^t \cdot b\f$, where  \f$a^t\f$ is the
   * transpose of \f$a\f$.
   */
  for (unsigned i = 0; i < p; i++) {
    atb[i] = new double[3];
    for (unsigned k = 0; k < 3; k++) {
      atb[i][k] = 0;
      for (unsigned j = 0; j < n; j++) {
        atb[i][k] += a[j][i] * b[(3 * j) + k];
      }
    }
  }

  return;
}

} // namespace pps
