/**
 * \file ppsfrompnt.cpp
 *
 * \brief  Implementation of  the class  PPSfromPNT that  represents a
 * Parametric Pseudo Surface, which approximates a PN triangle surface
 * defined on a  triangle mesh represented by a  Doubly Connected List
 * (DCEL).
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

#include "ppsfrompnt.hpp" // PPSfromPNT

#include <exception.hpp> // ERROR_UNLESS

#include <cmath> // std::fabs, std::sqrt

namespace ppsfrompnt {

using pps::Bezier;
using pps::PPS;

PPSfromPNT::PPSfromPNT(Surface *mesh) : PPS<Surface>(mesh) {
  //
  // Sets the owner of each half-edge attribute set.
  //
  for (EdgeIterator eit = edges_begin(); !is_done(eit); move_forward(eit)) {
    Edge *edge = get_edge(eit);

    Halfedge *h1 = edge->get_first_halfedge();
    Halfedge *h2 = edge->get_second_halfedge();

    h1->get_attributes().set_owner(h1);
    h2->get_attributes().set_owner(h2);
  }

  //
  // Create the PN triangle surface.
  //
  build_pnt_surface();
}

void PPSfromPNT::eval_surface(Face *face, double u, double v, double w,
                              double &x, double &y, double &z) const {
  //
  // Make sure the face has a PN triangle associated with it.
  //
  PNTriangle *patch = face->get_attributes().get_patch();

  ERROR_UNLESS(patch != 0, "Pointer to Bezier patch is null");

  //
  // Make sure the coordinates are non-negative and add up to 1.
  //
  ERROR_UNLESS((u >= 0) && (u <= 1) && (v >= 0) && (v <= 1) && (w >= 0) &&
                   (w <= 1) && (std::abs(1 - (u + v + w)) <= 1e-15),
               "Barycentric coordinates are invalid");
  //
  // Evaluate the patch.
  //
  patch->point(u, v, x, y, z);
}

bool PPSfromPNT::mesh_has_boundary() const { return false; }

bool PPSfromPNT::mesh_is_simplicial() const { return true; }

unsigned int PPSfromPNT::get_id(Halfedge *h) const {
  return h->get_attributes().get_pps_id();
}

typename PPSfromPNT::Vertex *PPSfromPNT::get_org(Halfedge *h) const {
  return h->get_origin();
}

typename PPSfromPNT::Vertex *PPSfromPNT::get_dst(Halfedge *h) const {
  return h->get_next()->get_origin();
}

typename PPSfromPNT::Edge *PPSfromPNT::get_edge(Halfedge *h) const {
  return h->get_edge();
}

typename PPSfromPNT::Face *PPSfromPNT::get_face(Halfedge *h) const {
  return h->get_face();
}

typename PPSfromPNT::Halfedge *PPSfromPNT::get_prev(Halfedge *h) const {
  return h->get_prev();
}

typename PPSfromPNT::Halfedge *PPSfromPNT::get_next(Halfedge *h) const {
  return h->get_next();
}

typename PPSfromPNT::Halfedge *PPSfromPNT::get_mate(Halfedge *h) const {
  return h->get_mate();
}

typename PPSfromPNT::Halfedge *PPSfromPNT::get_halfedge(Face *face) const {
  return face->get_halfedge();
}

typename PPSfromPNT::Halfedge *PPSfromPNT::get_halfedge(Vertex *vertex) const {
  return vertex->get_halfedge();
}

unsigned PPSfromPNT::get_degree(Vertex *vertex) const {
  return vertex->get_halfedge()->get_attributes().get_origin_vertex_degree();
}

Bezier *PPSfromPNT::get_shape_function(Vertex *vertex) const {
  return vertex->get_attributes().get_patch();
}

void PPSfromPNT::set_shape_function(Vertex *vertex, Bezier *patch) {
  vertex->get_attributes().set_patch(patch);
}

typename PPSfromPNT::VertexIterator PPSfromPNT::vertices_begin() const {
  return get_mesh()->vertices_begin();
}

bool PPSfromPNT::is_done(const VertexIterator &iterator) const {
  return iterator == get_mesh()->vertices_end();
}

void PPSfromPNT::move_forward(VertexIterator &iterator) const { ++iterator; }

typename PPSfromPNT::Vertex *
PPSfromPNT::get_vertex(const VertexIterator &iterator) const {
  return *iterator;
}

typename PPSfromPNT::EdgeIterator PPSfromPNT::edges_begin() const {
  return get_mesh()->edges_begin();
}

bool PPSfromPNT::is_done(const EdgeIterator &iterator) const {
  return iterator == get_mesh()->edges_end();
}

void PPSfromPNT::move_forward(EdgeIterator &iterator) const { ++iterator; }

typename PPSfromPNT::Edge *
PPSfromPNT::get_edge(const EdgeIterator &iterator) const {
  return *iterator;
}

typename PPSfromPNT::FaceIterator PPSfromPNT::faces_begin() const {
  return get_mesh()->faces_begin();
}

bool PPSfromPNT::is_done(const FaceIterator &iterator) const {
  return iterator == get_mesh()->faces_end();
}

void PPSfromPNT::move_forward(FaceIterator &iterator) const { ++iterator; }

typename PPSfromPNT::Face *
PPSfromPNT::get_face(const FaceIterator &iterator) const {
  return *iterator;
}

void PPSfromPNT::build_pnt_surface() {
  //
  // For each face of the given triangle mesh, compute a PN triangle.
  //
  for (FaceIterator fit = faces_begin(); !is_done(fit); move_forward(fit)) {
    Face *face = get_face(fit);

    Halfedge *h1 = face->get_halfedge();
    Halfedge *h2 = h1->get_next();
    Halfedge *h3 = h1->get_prev();

    Vertex *v1 = h1->get_origin();
    Vertex *v2 = h2->get_origin();
    Vertex *v3 = h3->get_origin();

    double p1[3];
    double p2[3];
    double p3[3];

    p1[0] = v1->x();
    p1[1] = v1->y();
    p1[2] = v1->z();

    p2[0] = v2->x();
    p2[1] = v2->y();
    p2[2] = v2->z();

    p3[0] = v3->x();
    p3[1] = v3->y();
    p3[2] = v3->z();

    double n1[3];
    double n2[3];
    double n3[3];

    compute_vertex_normal_vector(v1, n1[0], n1[1], n1[2]);
    compute_vertex_normal_vector(v2, n2[0], n2[1], n2[2]);
    compute_vertex_normal_vector(v3, n3[0], n3[1], n3[2]);

    //
    // Assign the PN triangle to the face
    //
    face->get_attributes().set_patch(new PNTriangle(p1, p2, p3, n1, n2, n3));
  }

  return;
}

void PPSfromPNT::compute_vertex_normal_vector(Vertex *vertex, double &x,
                                              double &y, double &z) {
  unsigned i = 0;

  x = 0;
  y = 0;
  z = 0;

  Halfedge *h1 = vertex->get_halfedge();
  Halfedge *h2 = h1;

  do {
    double xaux;
    double yaux;
    double zaux;

    compute_face_normal_vector(h2->get_face(), xaux, yaux, zaux);

    x += xaux;
    y += yaux;
    z += zaux;

    ++i;

    h2 = h2->get_prev()->get_mate();
  } while (h2 != h1);

  x /= double(i);
  y /= double(i);
  z /= double(i);

  double length = std::sqrt((x * x) + (y * y) + (z * z));

  x /= length;
  y /= length;
  z /= length;

  return;
}

void PPSfromPNT::compute_face_normal_vector(Face *face, double &x, double &y,
                                            double &z) {
  Halfedge *h1 = face->get_halfedge();
  Halfedge *h2 = h1->get_next();
  Halfedge *h3 = h1->get_prev();

  double x1 = h1->get_origin()->x();
  double y1 = h1->get_origin()->y();
  double z1 = h1->get_origin()->z();

  double x2 = h2->get_origin()->x();
  double y2 = h2->get_origin()->y();
  double z2 = h2->get_origin()->z();

  double x3 = h3->get_origin()->x();
  double y3 = h3->get_origin()->y();
  double z3 = h3->get_origin()->z();

  x2 -= x1;
  y2 -= y1;
  z2 -= z1;

  x3 -= x1;
  y3 -= y1;
  z3 -= z1;

  x = y2 * z3 - z2 * y3;
  y = z2 * x3 - x2 * z3;
  z = x2 * y3 - y2 * x3;

  double length = std::sqrt((x * x) + (y * y) + (z * z));

  x /= length;
  y /= length;
  z /= length;

  return;
}

} // namespace ppsfrompnt
