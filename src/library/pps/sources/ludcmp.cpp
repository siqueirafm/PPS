/**
 * \file ludcmp.cpp
 *
 * \brief  Implementation  of the  class  LUdcmp,  which represents  a
 * linear  system solver  based on  LU decomposition.   This  class is
 * basically a wrapper for the  LU decomposition code described in the
 * following book:
 *
 * William H. Press; Saul.  A. Teukolsky; William T. Vetterling; Brian
 * P.   Flannery.   Numerical Recipes  in  C:  The  Art of  Scientific
 * Computing. Second Edition, Cambridge University Press, 1992.
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

#include "ludcmp.hpp"

#include <exception.hpp> // ERROR_UNLESS

#include <cmath> // std::abs

namespace pps {

LUdcmp::LUdcmp(double **mat, unsigned n) {
  _nele = n;
  _matA = allocate(_nele);

  for (unsigned i = 0; i < _nele; i++) {
    for (unsigned j = 0; j < _nele; j++) {
      _matA[i][j] = mat[i][j];
    }
  }

  _perm.resize(_nele);
  _sign = 1;

  decomp();

  return;
}

LUdcmp::LUdcmp(const LUdcmp &lu) : _nele(lu._nele) {
  _matA = allocate(_nele);

  _perm.resize(_nele);
  _sign = lu._sign;

  for (unsigned i = 0; i < _nele; i++) {
    _perm[i] = lu._perm[i];

    for (unsigned j = 0; j < _nele; j++) {
      _matA[i][j] = lu._matA[i][j];
    }
  }

  return;
}

LUdcmp::~LUdcmp() {
  deallocate(_matA, _nele);

  return;
}

void LUdcmp::solve(double *b, unsigned n) const {
  ERROR_UNLESS(n == _nele, "Invalid number of array elements");

  int i, ii = 0, ip, j;
  double sum;
  for (i = 0; i < (int)n; i++) {
    ip = _perm[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii != 0) {
      for (j = ii - 1; j < i; j++) {
        sum -= _matA[i][j] * b[j];
      }
    } else if (sum != 0.0) {
      ii = i + 1;
    }

    b[i] = sum;
  }

  for (i = n - 1; i >= 0; i--) {
    sum = b[i];
    for (j = i + 1; j < (int)n; j++) {
      sum -= _matA[i][j] * b[j];
    }

    b[i] = sum / _matA[i][i];
  }

  return;
}

double **LUdcmp::allocate(unsigned n) const {
  double **mat = (double **)new double *[n];

  for (unsigned i = 0; i < n; i++) {
    mat[i] = (double *)new double[n];
  }

  return mat;
}

void LUdcmp::deallocate(double **mat, unsigned n) const {
  for (unsigned i = 0; i < n; i++) {
    delete[] mat[i];
  }

  delete mat;

  return;
}

void LUdcmp::decomp() {
  const double TINY = 1.0e-20;
  int i, imax = 0, j, k;
  double big, dum, sum, temp;

  int n = (int)_nele;
  std::vector<double> vv(n);
  _sign = 1;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++) {
      if ((temp = std::abs(_matA[i][j])) > big) {
        big = temp;
      }
    }

    ERROR_UNLESS(big != 0.0, "Failed to compute LU decomposition");

    vv[i] = 1.0 / big;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = _matA[i][j];

      for (k = 0; k < i; k++) {
        sum -= _matA[i][k] * _matA[k][j];
      }

      _matA[i][j] = sum;
    }

    big = 0.0;

    for (i = j; i < n; i++) {
      sum = _matA[i][j];
      for (k = 0; k < j; k++) {
        sum -= _matA[i][k] * _matA[k][j];
      }

      _matA[i][j] = sum;

      if ((dum = vv[i] * std::abs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }

    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = _matA[imax][k];
        _matA[imax][k] = _matA[j][k];
        _matA[j][k] = dum;
      }

      _sign = -_sign;
      vv[imax] = vv[j];
    }

    _perm[j] = imax;

    if (_matA[j][j] == 0.0) {
      _matA[j][j] = TINY;
    }

    if (j != (n - 1)) {
      dum = 1.0 / (_matA[j][j]);
      for (i = j + 1; i < n; i++) {
        _matA[i][j] *= dum;
      }
    }
  }

  return;
}

} // namespace pps
