/**
 * \file loop_eval.cpp
 *
 * \brief Implementation of class \c EvalStruct, which represents a Loop
 * subdivision surface exact evaluator.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 1.0
 * \date December 2008
 */

#include "loop_eval.hpp"

#include <cmath>   // pow, log2, floor
#include <cstdlib> // exit
#include <stdio.h> // fopen, fclose, fread, fwrite, fprintf

/*********************************************************************
 *
 * This  code  was borrowed  from  Jos  Stam,  who created  the  first
 * approach for exact evaluation of Catmull-Clark and Loop subdivision
 * surfaces:
 *
 * Stam, Jos. "Exact  Evaluation of Catmull-Clark Subdivision Surfaces
 * at Arbitrary Parameter  Values".  In Computer Graphics Proceedings,
 * Annual Conference Series, 1998, pages 395-404, July 1998.
 *
 *
 * The  class EvalStruct  is basically  a C++  wrapper class  for the
 * original C code by Stam, who kindly provided me with his original C
 * implementation.  If  you find any  problems with this  code, please
 * email ME  and not  him (as I  am the  one responsible for  the code
 * below):
 *
 * mfsiqueira at gmail dot com
 *
 ********************************************************************/

namespace ppsfromloop {

EvalStruct::EvalStruct(const char *fname) : BIG_ONE(1e10) {
  FILE *f_in;
  int i, N, K, nread, nread1, nread2, nread3;

  f_in = fopen(fname, "r");

  if (!f_in) {
    fprintf(stderr, "read_eval_struct: cannot open %s for read\n", fname);
    exit(0);
  }

  nread = fread(&Nmax, sizeof(int), 1, f_in);

  if (nread != 1 || Nmax < 3) {
    fprintf(stderr, "read_eval_struct: wrong file format\n");
    fclose(f_in);
    exit(0);
  }

  ev_struct = new EVALSTRUCT *[Nmax - 2];

  if (!ev_struct) {
    fprintf(stderr, "read_eval_struct: cannot allocate memory\n");
    fclose(f_in);
    exit(0);
  }

  for (i = 0; i < Nmax - 2; i++)
    ev_struct[i] = NULL;

  for (i = 0; i < Nmax - 2; i++) {
    N = i + 3;
    K = N + 6;
    ev_struct[i] = new EVALSTRUCT;
    if (!ev_struct[i]) {
      fprintf(stderr, "read_eval_struct: cannot allocate memory\n");
      fclose(f_in);
      exit(0);
    }

    ev_struct[i]->val = new double[K];
    ev_struct[i]->vecI = new double[K * K];
    ev_struct[i]->Phi = new double *[3];

    if (!ev_struct[i]->val || !ev_struct[i]->vecI || !ev_struct[i]->Phi) {
      fprintf(stderr, "read_eval_struct: cannot allocate memory\n");
      fclose(f_in);
      exit(0);
    }

    ev_struct[i]->Phi[0] = new double[K * 12];
    ev_struct[i]->Phi[1] = new double[K * 12];
    ev_struct[i]->Phi[2] = new double[K * 12];

    if (!ev_struct[i]->Phi[0] || !ev_struct[i]->Phi[1] ||
        !ev_struct[i]->Phi[2]) {
      fprintf(stderr, "read_eval_struct: cannot allocate memory\n");
      fclose(f_in);
      exit(0);
    }

    nread1 = fread(ev_struct[i]->val, sizeof(double), K, f_in);
    nread3 = fread(ev_struct[i]->vecI, sizeof(double), K * K, f_in);

    if (nread1 != K || nread3 != K * K) {
      fprintf(stderr, "read_eval_struct: error in data file %s\n", fname);
      fclose(f_in);
      exit(0);
    }

    nread1 = fread(ev_struct[i]->Phi[0], sizeof(double), K * 12, f_in);
    nread2 = fread(ev_struct[i]->Phi[1], sizeof(double), K * 12, f_in);
    nread3 = fread(ev_struct[i]->Phi[2], sizeof(double), K * 12, f_in);

    if (nread1 != K * 12 || nread2 != K * 12 || nread3 != K * 12) {
      fprintf(stderr, "read_eval_struct: error in data file %s\n", fname);
      fclose(f_in);
      exit(0);
    }
  }

  fclose(f_in);
}

EvalStruct::~EvalStruct() {
  int i;

  if (!ev_struct)
    return;

  for (i = 0; i < Nmax - 2; i++) {
    if (!ev_struct[i])
      continue;

    if (ev_struct[i]->val)
      delete[] ev_struct[i]->val;

    if (ev_struct[i]->vecI)
      delete[] ev_struct[i]->vecI;

    if (ev_struct[i]->Phi) {
      if (ev_struct[i]->Phi[0])
        delete[] ev_struct[i]->Phi[0];

      if (ev_struct[i]->Phi[1])
        delete[] ev_struct[i]->Phi[1];

      if (ev_struct[i]->Phi[2])
        delete[] ev_struct[i]->Phi[2];
    }

    delete ev_struct[i];
  }

  delete[] ev_struct;
}

void EvalStruct::eval_surf(EVALTYPE der, double *pC, int N, double v, double w,
                           double s[3][3]) const {
  EVALSTRUCT *ev;
  int i, j, m, k, n, p, K, p2;
  double *Phi;
  double bs, pow_l, pow_l0 = 0;
  double v0, w0;
  double b[3][12];

  if (N < 3) {
    fprintf(stderr, "eval_surf: N should be no less than 3\n");
    return;
  }

  if (v < 0.0001 && w < 0.0001)
    v = w = 0.0001;

  K = N + 6;
  ev = ev_struct[N - 3];

  n = (int)std::floor(1 - std::log2(v + w));

  p2 = (n == 0) ? 1 : 1 << (n - 1);

  v0 = v * p2;
  w0 = w * p2;

  if (v0 > 0.5) {
    v0 = 2 * v0 - 1;
    w0 = 2 * w0;
    Phi = ev->Phi[0];
    k = 0;
  } else if (w0 > 0.5) {
    v0 = 2 * v0;
    w0 = 2 * w0 - 1;
    Phi = ev->Phi[2];
    k = 2;
  } else {
    v0 = 1 - 2 * v0;
    w0 = 1 - 2 * w0;
    Phi = ev->Phi[1];
    k = 1;
  }

  switch (der) {
  case EVAL_VALUE:
    eval_T_spline(b[0], v0, w0);
    for (p = 0; p < 3; p++) {
      s[0][p] = 0.0;
    }
    break;

  case EVAL_DER:
    eval_T_spline_v(b[0], v0, w0);
    eval_T_spline_w(b[1], v0, w0);
    for (p = 0; p < 3; p++) {
      s[0][p] = s[1][p] = 0.0;
    }
    break;

  case EVAL_DER2:
    eval_T_spline_vv(b[0], v0, w0);
    eval_T_spline_vw(b[1], v0, w0);
    eval_T_spline_ww(b[2], v0, w0);
    for (p = 0; p < 3; p++) {
      s[0][p] = s[1][p] = s[2][p] = 0.0;
    }
    break;
  }

  for (i = 0; i < K; i++) {
    for (m = 0; m < der; m++) {
      bs = 0.0;
      for (j = 0; j < 12; j++) {
        bs += Phi[IX(i, j, K)] * b[m][j];
      }

      pow_l = n <= 1 ? 1.0 : std::pow(ev->val[i], (double)n - 1);

      if (N == 3 && i == K - 2) {
        pow_l0 = n <= 1 ? 0.0 : n == 2 ? 1.0 : (n - 2) * pow_l / ev->val[i];
      }

      for (p = 0; p < 3; p++) {
        s[m][p] += pow_l * pC[IX(i, p, K)] * bs;

        if (N == 3 && i == K - 2) {
          s[m][p] += pow_l0 * pC[IX(i + 1, p, K)] * bs;
        }
      }
    }
  }

  switch (der) {
  case EVAL_VALUE:
    break;

  case EVAL_DER:
    for (p = 0; p < 3; p++) {
      s[0][p] *= k == 1 ? -2 * p2 : 2 * p2;
      s[1][p] *= k == 1 ? -2 * p2 : 2 * p2;
    }
    break;

  case EVAL_DER2:
    for (p = 0; p < 3; p++) {
      s[0][p] *= 4 * p2 * p2;
      s[1][p] *= 4 * p2 * p2;
      s[2][p] *= 4 * p2 * p2;
    }
    break;
  }
}

double *EvalStruct::proj_control_points(double *C, int N) const {
  EVALSTRUCT *ev;
  double *pC, sum;
  int i, j, k, K;

  if (N < 3) {
    fprintf(stderr, "proj_control_points: N should be bigger than 3\n");
    return (NULL);
  }

  ev = ev_struct[N - 3];

  K = N + 6;

  pC = get_mem(K * 3);
  if (!pC) {
    fprintf(stderr, "proj_control_points: cannot allocate memory\n");
    return (NULL);
  }

  for (k = 0; k < 3; k++) {

    for (i = 0; i < K; i++) {
      sum = 0.0;
      for (j = 0; j < K; j++) {
        sum += ev->vecI[IX(i, j, K)] * C[IX(j, k, K)];
      }

      pC[IX(i, k, K)] = sum;
    }
  }

  return (pC);
}

void EvalStruct::printout() const {
  int i, j, k, l, N, K;

  printf("Nmax = %d\n\n", Nmax);

  for (i = 0; i < Nmax - 2; i++) {
    N = i + 3;
    K = N + 6;

    printf("N=%d\n\n", N);

    printf("Eigenvalues:\n\n");

    for (j = 0; j < K; j++)
      printf("%6.3f ", ev_struct[i]->val[j]);

    printf("\n\n");

    printf("Inverse of Eigenvectors:\n\n");

    for (j = 0; j < K; j++) {
      for (k = 0; k < K; k++)
        printf("%6.3f ", ev_struct[i]->vecI[IX(j, k, K)]);

      printf("\n");
    }

    printf("\n\n");

    printf("Coefficients of the Eigenbasis functions\n\n");

    for (k = 0; k < 3; k++) {
      printf("k=%d :\n\n", k);

      for (j = 0; j < K; j++) {
        for (l = 0; l < 12; l++)
          printf("%6.3f ", ev_struct[i]->Phi[k][IX(j, l, K)]);

        printf("\n");
      }

      printf("\n\n");
    }

    printf("\n\n");
  }
}

int EvalStruct::get_max_degree() const { return Nmax; }

int EvalStruct::IX(int i, int j, int n) const { return i + n * j; }

double *EvalStruct::get_mem(int size) const {
  double *d;

  d = new double[size];

  return (d);
}

void EvalStruct::eval_T_spline(double *b, double v, double w) const {
  double u;
  double t1, t2, t3, t4, t6, t9, t10, t11, t13, t14, t15, t16, t17, t19, t20,
      t21, t22, t23, t25, t26, t27, t28, t31, t32;

  u = 1 - v - w;

  t1 = u * u;
  t2 = t1 * t1;
  t3 = t1 * u;
  t4 = t3 * v;
  t6 = t3 * w;
  t9 = t1 * v * w;
  t10 = v * v;
  t11 = t1 * t10;
  t13 = u * t10 * w;
  t14 = t10 * v;
  t15 = u * t14;
  t16 = t14 * w;
  t17 = t10 * t10;
  t19 = w * w;
  t20 = t1 * t19;
  t21 = t19 * w;
  t22 = u * t21;
  t23 = t19 * t19;
  t25 = u * v * t19;
  t26 = v * t21;
  t27 = t10 * t19;
  t28 = t2 / 2 + 2.0 * t6 + 2.0 * t20 + 2.0 / 3.0 * t22 + t23 / 12 + 2.0 * t4 +
        5.0 * t9 + 3.0 * t25 + t26 / 2 + 2.0 * t11 + 3.0 * t13 + t27 +
        2.0 / 3.0 * t15 + t16 / 2 + t17 / 12;
  t31 = t2 / 12 + t6 / 2 + t20 + t22 / 2 + t23 / 12 + 2.0 / 3.0 * t4 +
        3.0 * t9 + 3.0 * t25 + 2.0 / 3.0 * t26 + 2.0 * t11 + 5.0 * t13 +
        2.0 * t27 + 2.0 * t15 + 2.0 * t16 + t17 / 2;
  t32 = t2 / 12 + 2.0 / 3.0 * t6 + 2.0 * t20 + 2.0 * t22 + t23 / 2 + t4 / 2 +
        3.0 * t9 + 5.0 * t25 + 2.0 * t26 + t11 + 3.0 * t13 + 2.0 * t27 +
        t15 / 2 + 2.0 / 3.0 * t16 + t17 / 12;

  b[0] = t2 / 12 + t4 / 6;
  b[1] = t2 / 12 + t6 / 6;
  b[2] = t2 / 12 + t6 / 6 + t4 / 2 + t9 / 2 + t11 + t13 / 2 + t15 / 2 +
         t16 / 6 + t17 / 12;
  b[3] = t28;
  b[4] = t2 / 12 + t6 / 2 + t20 + t22 / 2 + t23 / 12 + t4 / 6 + t9 / 2 +
         t25 / 2 + t26 / 6;
  b[5] = t15 / 6 + t17 / 12;
  b[6] = t31;
  b[7] = t32;
  b[8] = t22 / 6 + t23 / 12;
  b[9] = t16 / 6 + t17 / 12;
  b[10] = t22 / 6 + t23 / 12 + t25 / 2 + t26 / 2 + t13 / 2 + t27 + t15 / 6 +
          t16 / 2 + t17 / 12;
  b[11] = t23 / 12 + t26 / 6;
}

void EvalStruct::eval_T_spline_v(double *bv, double v, double w) const {
  double u;
  double t1, t2, t3, t4, t5, t7, t8, t11, t12, t13, t14, t15, t17, t18, t19,
      t20, t21, t23, t24, t25, t26, t27, t29, t32;

  u = 1 - v - w;

  t1 = u * u;
  t2 = t1 * u;
  t3 = t2 / 6.0;
  t4 = t1 * v;
  t5 = t4 / 2.0;
  t7 = t2 / 3.0;
  t8 = t1 * w;
  t11 = v * v;
  t12 = u * t11;
  t13 = t12 / 2.0;
  t14 = t11 * v;
  t15 = t14 / 6.0;
  t17 = w * w;
  t18 = u * t17;
  t19 = t17 * w;
  t20 = t19 / 6.0;
  t21 = 2.0 * t4;
  t23 = u * v * w;
  t24 = 4.0 * t23;
  t25 = v * t17;
  t26 = 2.0 * t12;
  t27 = t11 * w;
  t29 = t14 / 3.0;
  t32 = t19 / 3.0;

  bv[0] = -t3 - t5;
  bv[1] = -t7 - t8 / 2.0;
  bv[2] = t3 + t5 - t13 - t15;
  bv[3] = -t8 - t18 - t20 - t21 - t24 - t25 - t26 - 3.0 / 2.0 * t27 - t29;
  bv[4] = -t3 - t8 - 3.0 / 2.0 * t18 - t32 - t5 - t23 - t25 / 2.0;
  bv[5] = t15 + t13;
  bv[6] = t7 + 3.0 / 2.0 * t8 + t18 + t20 + t21 + t24 + t25 + t26 + t27;
  bv[7] = t3 + t8 + t18 + t5 - t25 - t13 - t27 - t15;
  bv[8] = -t20;
  bv[9] = t27 / 2.0 + t29;
  bv[10] = t32 + 3.0 / 2.0 * t25 + t18 / 2.0 + t27 + t23 + t15 + t13;
  bv[11] = t20;
}

void EvalStruct::eval_T_spline_w(double *bw, double v, double w) const {
  double u;
  double t1, t2, t3, t4, t7, t8, t9, t12, t13, t14, t16, t18, t19, t21, t22,
      t23, t24, t25, t26, t27, t28, t30, t32, t33;

  u = 1 - v - w;

  t1 = u * u;
  t2 = t1 * u;
  t3 = t2 / 3.0;
  t4 = t1 * v;
  t7 = t2 / 6.0;
  t8 = t1 * w;
  t9 = t8 / 2.0;
  t12 = u * v * w;
  t13 = v * v;
  t14 = u * t13;
  t16 = t13 * w;
  t18 = t13 * v;
  t19 = t18 / 3.0;
  t21 = 2.0 * t8;
  t22 = w * w;
  t23 = u * t22;
  t24 = 2.0 * t23;
  t25 = t22 * w;
  t26 = t25 / 3.0;
  t27 = 4.0 * t12;
  t28 = v * t22;
  t30 = t18 / 6.0;
  t32 = t23 / 2.0;
  t33 = t25 / 6.0;

  bw[0] = -t3 - t4 / 2.0;
  bw[1] = -t7 - t9;
  bw[2] = -t7 - t9 - t4 - t12 - 3.0 / 2.0 * t14 - t16 / 2.0 - t19;
  bw[3] = -t21 - t24 - t26 - t4 - t27 - 3.0 / 2.0 * t28 - t14 - t16 - t30;
  bw[4] = t7 + t9 - t32 - t33;
  bw[5] = -t30;
  bw[6] = t7 + t9 - t32 - t33 + t4 - t28 + t14 - t16;
  bw[7] = t3 + t21 + t24 + 3.0 / 2.0 * t4 + t27 + t28 + t14 + t16 + t30;
  bw[8] = t33 + t32;
  bw[9] = t30;
  bw[10] = t33 + t32 + t28 + t12 + 3.0 / 2.0 * t16 + t14 / 2.0 + t19;
  bw[11] = t26 + t28 / 2.0;
}

void EvalStruct::eval_T_spline_vv(double *bvv, double v, double w) const {
  double u;
  double t1, t2, t3, t8, t9, t11, t12;

  u = 1 - v - w;

  t1 = u * v;
  t2 = u * u;
  t3 = u * w;
  t8 = v * w;
  t9 = v * v;
  t11 = w * w;
  t12 = t3 + t11 + t1 + t8;

  bvv[0] = t1;
  bvv[1] = t2 + t3;
  bvv[2] = -2.0 * t1;
  bvv[3] = -2.0 * t2 - 2.0 * t3 + t8 + t9;
  bvv[4] = t12;
  bvv[5] = t1;
  bvv[6] = t2 + t3 - 2.0 * t8 - 2.0 * t9;
  bvv[7] = -2.0 * t12;
  bvv[8] = 0.0;
  bvv[9] = t8 + t9;
  bvv[10] = t12;
  bvv[11] = 0.0;
}

void EvalStruct::eval_T_spline_vw(double *bvw, double v, double w) const {
  double u;
  double t1, t2, t3, t5, t7, t8, t10, t11, t13;

  u = 1 - v - w;

  t1 = u * u;
  t2 = t1 / 2.0;
  t3 = u * v;
  t5 = u * w;
  t7 = v * v;
  t8 = t7 / 2.0;
  t10 = w * w;
  t11 = t10 / 2.0;
  t13 = 2.0 * v * w;

  bvw[0] = t2 + t3;
  bvw[1] = t2 + t5;
  bvw[2] = -t2 - t3 + t8;
  bvw[3] = t8 - t1 + t11 + t13;
  bvw[4] = -t2 - t5 + t11;
  bvw[5] = -t8;
  bvw[6] = -t7 + t2 - t11 - t13 - t5;
  bvw[7] = -t8 + t2 - t10 - t13 - t3;
  bvw[8] = -t11;
  bvw[9] = t8;
  bvw[10] = t11 + t13 + t5 + t8 + t3;
  bvw[11] = t11;
}

void EvalStruct::eval_T_spline_ww(double *bww, double v, double w) const {
  double u;
  double t1, t2, t4, t5, t6, t7, t9;

  u = 1 - v - w;

  t1 = u * u;
  t2 = u * v;
  t4 = u * w;
  t5 = v * w;
  t6 = v * v;
  t7 = t4 + t2 + t5 + t6;
  t9 = w * w;

  bww[0] = t1 + t2;
  bww[1] = t4;
  bww[2] = t7;
  bww[3] = -2.0 * t1 + t9 - 2.0 * t2 + t5;
  bww[4] = -2.0 * t4;
  bww[5] = 0.0;
  bww[6] = -2.0 * t7;
  bww[7] = t1 - 2.0 * t9 + t2 - 2.0 * t5;
  bww[8] = t4;
  bww[9] = 0.0;
  bww[10] = t7;
  bww[11] = t9 + t5;
}

} // namespace ppsfromloop
