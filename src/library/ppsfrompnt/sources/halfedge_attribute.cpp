/**
 * \file halfedge_attribute.cpp
 *
 * \brief Implementation of the class \c HalfedgeAttribute.
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

#include "halfedge_attribute.hpp" // HalfedgeAttribute

#include <exception.hpp> // ERROR_UNLESS

namespace ppsfrompnt {

HalfedgeAttribute::HalfedgeAttribute() {
  set_owner(0);
  set_pps_id(0);
  set_pps_id_flag(false);
  set_origin_vertex_degree(0);
  set_degree_flag(false);
}

HalfedgeAttribute::HalfedgeAttribute(Halfedge *h) {
  set_owner(h);
  set_pps_id(0);
  set_pps_id_flag(false);
  set_origin_vertex_degree(0);
  set_degree_flag(false);
}

typename HalfedgeAttribute::Halfedge *HalfedgeAttribute::get_owner() const {
  return _owner;
}

void HalfedgeAttribute::set_owner(Halfedge *h) { _owner = h; }

unsigned HalfedgeAttribute::get_pps_id() {
  if (get_pps_id_flag()) {
    return _id;
  }

  /*
   * To speed up  computation, compute the ID only  once and store
   * it  in  a  data  member.   Since  the  PPS  is  static,  this
   * optimization is safe.
   */

  set_pps_id(compute_pps_id());

  set_pps_id_flag(true);

  return _id;
}

void HalfedgeAttribute::set_pps_id(unsigned id) { _id = id; }

bool HalfedgeAttribute::get_pps_id_flag() const { return _id_flag; }

void HalfedgeAttribute::set_pps_id_flag(bool flag) { _id_flag = flag; }

unsigned HalfedgeAttribute::get_origin_vertex_degree() {
  if (get_degree_flag()) {
    return _degree;
  }

  /*
   * To speed up  computation, compute the degree only  once and store
   * it in a data member.  Since  the PPS is static, this optimization
   * is safe.
   */

  set_origin_vertex_degree(compute_origin_vertex_degree());

  set_degree_flag(true);

  return _degree;
}

void HalfedgeAttribute::set_origin_vertex_degree(unsigned degree) {
  _degree = degree;
}

bool HalfedgeAttribute::get_degree_flag() const { return _degree_flag; }

void HalfedgeAttribute::set_degree_flag(bool flag) { _degree_flag = flag; }

unsigned HalfedgeAttribute::compute_pps_id() const {
  ERROR_UNLESS(get_owner() != 0, "Pointer to attribute owner is null");

  unsigned i = 0;
  Halfedge *h = get_owner()->get_origin()->get_halfedge();

  while (h != get_owner()) {
    i++;
    h = h->get_prev()->get_mate();
  }

  return i;
}

unsigned HalfedgeAttribute::compute_origin_vertex_degree() const {
  ERROR_UNLESS(get_owner() != 0, "Pointer to attribute owner is null");

  unsigned i = 0;
  Halfedge *h = get_owner();
  do {
    ++i;

    h = h->get_prev()->get_mate();
  } while (h != get_owner());

  return i;
}

} // namespace ppsfrompnt
