/**
 * \file ppsfromloop.cpp
 *
 * \brief  Definition  of  the  class PPSfromLOOP  that  represents  a
 * Parametric  Pseudo Surface, which  approximates a  Loop subdivision
 * surface  defined  on  a  triangle  mesh  represented  by  a  Doubly
 * Connected List (DCEL).
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

#include "ppsfromloop.hpp" // PPSfromLOOP

#include <exception.hpp> // ERROR_UNLESS

#include <cmath> // std::abs

namespace ppsfromloop {

using pps::PPS;

namespace fs = std::filesystem;

PPSfromLOOP::PPSfromLOOP(Surface *mesh, const std::filesystem::path &filepath)
    : PPS<Surface>(mesh) {
  //
  // Create the Loop surface evaluator.
  //

  //
  // Load in the evaluator data table.
  //
  ERROR_UNLESS(fs::exists(filepath),
               "Given Loop surface evaluator table file cannot be found: " +
                   filepath.string());
  ERROR_UNLESS(
      fs::is_regular_file(filepath),
      "Given Loop surface evaluator table file is not a regular file: " +
          filepath.string());

  _evstruct = new EvalStruct(filepath.string().c_str());

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
  // Create the patches associated with the mesh base faces.
  //
  build_surface_evaluator();

  return;
}

PPSfromLOOP::~PPSfromLOOP() {
  if (_evstruct != 0) {
    delete _evstruct;
  }

  return;
}

void PPSfromLOOP::eval_surface(Face *face, double u, double v, double w,
                               double &x, double &y, double &z) const {
  //
  // Make  sure  the  face  has  a Loop  subdivision  surface  patch
  // associated with it.
  //
  LOOPatch *patch = face->get_attributes().get_patch();

  ERROR_UNLESS(patch, "Pointer to surface patch is the null pointer");

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

  return;
}

bool PPSfromLOOP::mesh_has_boundary() const { return false; }

bool PPSfromLOOP::mesh_is_simplicial() const { return true; }

unsigned int PPSfromLOOP::get_id(Halfedge *h) const {
  return h->get_attributes().get_pps_id();
}

typename PPSfromLOOP::Vertex *PPSfromLOOP::get_org(Halfedge *h) const {
  return h->get_origin();
}

typename PPSfromLOOP::Vertex *PPSfromLOOP::get_dst(Halfedge *h) const {
  return h->get_next()->get_origin();
}

typename PPSfromLOOP::Edge *PPSfromLOOP::get_edge(Halfedge *h) const {
  return h->get_edge();
}

typename PPSfromLOOP::Face *PPSfromLOOP::get_face(Halfedge *h) const {
  return h->get_face();
}

typename PPSfromLOOP::Halfedge *PPSfromLOOP::get_prev(Halfedge *h) const {
  return h->get_prev();
}

typename PPSfromLOOP::Halfedge *PPSfromLOOP::get_next(Halfedge *h) const {
  return h->get_next();
}

typename PPSfromLOOP::Halfedge *PPSfromLOOP::get_mate(Halfedge *h) const {
  return h->get_mate();
}

typename PPSfromLOOP::Halfedge *PPSfromLOOP::get_halfedge(Face *face) const {
  return face->get_halfedge();
}

typename PPSfromLOOP::Halfedge *
PPSfromLOOP::get_halfedge(Vertex *vertex) const {
  return vertex->get_halfedge();
}

unsigned PPSfromLOOP::get_degree(Vertex *vertex) const {
  return vertex->get_halfedge()->get_attributes().get_origin_vertex_degree();
}

Bezier *PPSfromLOOP::get_shape_function(Vertex *vertex) const {
  return vertex->get_attributes().get_patch();
}

void PPSfromLOOP::set_shape_function(Vertex *vertex, Bezier *patch) {
  vertex->get_attributes().set_patch(patch);
}

typename PPSfromLOOP::VertexIterator PPSfromLOOP::vertices_begin() const {
  return get_mesh()->vertices_begin();
}

bool PPSfromLOOP::is_done(const VertexIterator &iterator) const {
  return iterator == get_mesh()->vertices_end();
}

void PPSfromLOOP::move_forward(VertexIterator &iterator) const { ++iterator; }

typename PPSfromLOOP::Vertex *
PPSfromLOOP::get_vertex(const VertexIterator &iterator) const {
  return *iterator;
}

typename PPSfromLOOP::EdgeIterator PPSfromLOOP::edges_begin() const {
  return get_mesh()->edges_begin();
}

bool PPSfromLOOP::is_done(const EdgeIterator &iterator) const {
  return iterator == get_mesh()->edges_end();
}

void PPSfromLOOP::move_forward(EdgeIterator &iterator) const { ++iterator; }

typename PPSfromLOOP::Edge *
PPSfromLOOP::get_edge(const EdgeIterator &iterator) const {
  return *iterator;
}

typename PPSfromLOOP::FaceIterator PPSfromLOOP::faces_begin() const {
  return get_mesh()->faces_begin();
}

bool PPSfromLOOP::is_done(const FaceIterator &iterator) const {
  return iterator == get_mesh()->faces_end();
}

void PPSfromLOOP::move_forward(FaceIterator &iterator) const { ++iterator; }

typename PPSfromLOOP::Face *
PPSfromLOOP::get_face(const FaceIterator &iterator) const {
  return *iterator;
}

void PPSfromLOOP::build_surface_evaluator() {
  /**
   * For each face  of the given triangle surface  mesh, compute the
   * surface patch associated with  it in the given Loop subdivision
   * surface.
   */
  for (FaceIterator fit = faces_begin(); !is_done(fit); move_forward(fit)) {

    // -------------------------------------------------------------
    // Find the extraordinary vertex of the current face, if any.
    //
    // Recall that there can be AT MOST ONE extraordinary vertex per
    // face.
    // -------------------------------------------------------------

    Face *face = get_face(fit);

    Halfedge *he1 = face->get_halfedge();
    Halfedge *he2 = he1->get_next();
    Halfedge *he3 = he1->get_prev();

    unsigned nv1 = he1->get_attributes().get_origin_vertex_degree();
    unsigned nv2 = he2->get_attributes().get_origin_vertex_degree();
    unsigned nv3 = he3->get_attributes().get_origin_vertex_degree();

    Halfedge *he;
    unsigned nv;
    LOOPatch::BCORDER bo;

    ERROR_UNLESS(((nv2 == 6) && (nv3 == 6)) || ((nv1 == 6) && (nv3 == 6)) ||
                     ((nv1 == 6) && (nv2 == 6)),
                 "Neighboring vertices have incompatible vertex degree");

    if ((nv2 == 6) && (nv3 == 6)) {
      nv = nv1;
      he = he1;
      bo = LOOPatch::UVW;
    } else if ((nv1 == 6) && (nv3 == 6)) {
      nv = nv2;
      he = he2;
      bo = LOOPatch::VWU;
    } else { // ( ( nv1 == 6 ) && ( nv2 == 6 ) )
      nv = nv3;
      he = he3;
      bo = LOOPatch::WUV;
    }

    const unsigned int MAX_DEGREE = unsigned(_evstruct->get_max_degree());

    ERROR_UNLESS(nv <= MAX_DEGREE, "Vertex degree violates degree upper bound");

    // Compute  the "projected  points" of  the  subdivision surface
    // patch  associated with  the current  face of  the  given base
    // mesh.

    const unsigned int K = nv + 6;

    Halfedge **vhe = new Halfedge *[K];

    vhe[0] = he;
    vhe[1] = he->get_next();
    he = he->get_mate()->get_next()->get_mate();

    for (unsigned int i = 2; i <= nv; i++) {
      vhe[i] = he;
      he = he->get_next()->get_mate();
    }

    he = vhe[1]->get_mate()->get_next()->get_mate();

    for (unsigned int i = nv + 1; i <= nv + 3; i++) {
      vhe[i] = he;
      he = he->get_next()->get_mate();
    }

    he = vhe[1]->get_mate()->get_prev()->get_mate()->get_prev();

    for (unsigned int i = nv + 4; i <= nv + 5; i++) {
      vhe[i] = he;
      he = he->get_mate()->get_prev();
    }

    double *control_pts = new double[K * 3];

    for (unsigned int i = 0; i < K; i++) {
      control_pts[LOOPatch::index(0, i, K)] = vhe[i]->get_origin()->x();
      control_pts[LOOPatch::index(1, i, K)] = vhe[i]->get_origin()->y();
      control_pts[LOOPatch::index(2, i, K)] = vhe[i]->get_origin()->z();
    }

    delete[] vhe;

    double *proj_pts = _evstruct->proj_control_points(control_pts, int(nv));

    delete[] control_pts;

    //
    // Assign the  Loop subdivision  surface patch with  the current
    // face.
    //
    face->get_attributes().set_patch(new LOOPatch(_evstruct, nv, bo, proj_pts));

    delete[] proj_pts;
  }

  return;
}

} // namespace ppsfromloop
