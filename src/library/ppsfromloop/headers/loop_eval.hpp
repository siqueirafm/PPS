#pragma once

/**
 * \file loop_eval.hpp
 *
 * \brief Definition of  class \c EvalStruct, which  represents a Loop
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
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

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
 * The  class EvalStruct  is basically  a  C++ wrapper  class for  the
 * original C code by Stam, who kindly provided me with his original C
 * implementation.  If  you find any  problems with this  code, please
 * email ME  and not  him (as I  am the  one responsible for  the code
 * below):
 *
 * mfsiqueira at gmail dot com
 *
 ********************************************************************/

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
 * \class EvalStruct
 *
 * \brief  This class  represents a  Loop subdivision  surface exact
 * evaluator.
 */
class EvalStruct {
public:
  enum EVALTYPE { EVAL_VALUE = 1, EVAL_DER = 2, EVAL_DER2 = 3 };

  struct EVALSTRUCT {
    double *val;
    double *vecI;
    double **Phi;
  };

  explicit EvalStruct(const char *);

  ~EvalStruct();

  double *proj_control_points(double *, int) const;

  void eval_surf(EVALTYPE, double *, int, double, double, double[3][3]) const;

  void printout() const;

  int get_max_degree() const;

private:
  int IX(int i, int j, int n) const;

  double *get_mem(int) const;

  void eval_T_spline(double *, double, double) const;
  void eval_T_spline_v(double *, double, double) const;
  void eval_T_spline_w(double *, double, double) const;
  void eval_T_spline_vv(double *, double, double) const;
  void eval_T_spline_vw(double *, double, double) const;
  void eval_T_spline_ww(double *, double, double) const;

private:
  const double BIG_ONE;
  EVALSTRUCT **ev_struct;
  int Nmax;
};

} // namespace ppsfromloop

/** @} */ // end of group class.
