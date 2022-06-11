#pragma once

/**
 * \file ppsfromloop.h
 *
 * \brief  Definition  of  class  \c  PPSfromLOOP  that  represents  a
 * Parametric Pseudo  Surface, which  approximates a  Loop subdivision
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
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#include "face_attribute.hpp"     // FaceAttribute
#include "halfedge_attribute.hpp" // HalfedgeAttribute
#include "loop_eval.hpp"          // EvalStruct
#include "loop_patch.hpp"         // LOOPatch
#include "vertex_attribute.hpp"   // VertexAttribute

#include <bezier.hpp>  // pps::Bezier
#include <pps.hpp>     // pps::PPS
#include <surface.hpp> // dcel::Surface

#include <filesystem> // std::filesystem::path

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

using pps::Bezier;
using pps::PPS;

/**
 * \class PPSfromLOOP
 *
 * \brief  This class  represents a  Parametric  Pseudo-Surface that
 * approximates a Loop subdivision surface.
 *
 */
class PPSfromLOOP : public PPS<dcel::Surface<VertexAttribute, FaceAttribute,
                                             int, HalfedgeAttribute>> {
public:
  /**
   * \using Halfedge
   *
   * \brief   Defines  Surface  as   an  alias   for  dcel::Surface<
   * VertexAttribute , FaceAttribute , int , HalfedgeAttribute >
   */
  using Surface =
      dcel::Surface<VertexAttribute, FaceAttribute, int, HalfedgeAttribute>;

  /**
   * \fn PPSfromLOOP( Surface* mesh , const std::filesystem::path& filepath )
   *
   * \brief Creates an instance of this PPSfromLOOP class.
   *
   * \param mesh A Pointer to a surface mesh represented by a DCEL.
   * \param filepath Full path to the file containing the Loop surface
   * evaluator data table.
   *
   */
  PPSfromLOOP(Surface *mesh, const std::filesystem::path &filepath);

  /**
   * \fn ~PPSfromLOOP()
   *
   * \brief Destroys an instance of this class.
   */
  ~PPSfromLOOP();

  /**
   * \fn void eval_surface( Face* face , double u , double v , double w ,
   * double& x , double& y , double& z ) const
   *
   * \brief Computes a point on the Loop's subdivision surface patch
   * associated  with a given  face of  the PPS  underlying triangle
   * mesh.
   *
   * \param face Pointer to one face of the underlying mesh.
   * \param u First barycentric coordinate of a point on the face.
   * \param v Second barycentric coordinate of a point on the face.
   * \param w Third barycentric coordinate of a point on the face.
   * \param x First Cartesian coordinate of a point on the PPS image.
   * \param y Second Cartesian coordinate of a point on the PPS image.
   * \param z Third Cartesian coordinate of a point on the PPS image.
   */
  void eval_surface(Face *face, double u, double v, double w, double &x,
                    double &y, double &z) const;

  /**
   * \fn bool mesh_has_boundary() const
   *
   * \brief  Determines if  the underlying  mesh of  this PPS  has a
   * non-empty boundary.
   *
   * \returns  The logic  value true  if  the mesh  had a  non-empty
   * boundary and the logic value false otherwise.
   */
  bool mesh_has_boundary() const;

  /**
   * \fn bool mesh_is_simplicial() const
   *
   * \brief  Determines if  the underlying  mesh  of this  PPS is  a
   * simplicial complex.
   *
   * \returns  The logic  value true  if  the mesh  is a  simplicial
   * complex and the logic value false otherwise.
   */
  bool mesh_is_simplicial() const;

  /**
   * \fn unsigned int get_id( Halfedge* h ) const
   *
   * \brief Returns the identifier of a given half-edge.
   *
   * \param h Pointer to a half-edge.
   *
   * \returns The identifier of the given half-edge.
   */
  unsigned int get_id(Halfedge *h) const;

  /**
   * \fn Vertex* get_org( Halfedge* h ) const
   *
   * \brief Returns the  origin vertex of a given  half-edge of this
   * PPS underlying mesh.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the origin vertex of a given half-edge.
   */
  Vertex *get_org(Halfedge *h) const;

  /**
   * \fn inline Vertex* get_dst( Halfedge* h ) const
   *
   * \brief Returns  the destination vertex of a  given half-edge of
   * this PPS underlying mesh.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns  A  pointer  to  the  destination vertex  of  a  given
   * half-edge.
   */
  Vertex *get_dst(Halfedge *h) const;

  /**
   * \fn inline Edge* get_edge( Halfedge* h ) const
   *
   * \brief  Returns  the  edge   a  given  half-edge  of  this  PPS
   * underlying mesh belongs to.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the edge a given half-edge belongs to.
   */
  Edge *get_edge(Halfedge *h) const;

  /**
   * \fn Face* get_face( Halfedge* h ) const
   *
   * \brief  Returns  the  face   a  given  half-edge  of  this  PPS
   * underlying mesh belongs to.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the face a given half-edge belongs to.
   */
  Face *get_face(Halfedge *h) const;

  /**
   * \fn Halfedge* get_prev( Halfedge* h ) const
   *
   * \brief Returns the half-edge that precedes a given half-edge of
   * this PPS underlying  mesh in the face cycle  of half-edges that
   * contains both half-edges.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns  A pointer  to  the half-edge  that  precedes a  given
   * half-edge in their common face half-edge cycle.
   */
  Halfedge *get_prev(Halfedge *h) const;

  /**
   * \fn Halfedge* get_next( Halfedge* h ) const
   *
   * \brief Returns the half-edge that succeeds a given half-edge of
   * this PPS underlying  mesh in the face cycle  of half-edges that
   * contains both half-edges.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns  A pointer  to  the half-edge  that  succeeds a  given
   * half-edge in their common face half-edge cycle.
   */
  Halfedge *get_next(Halfedge *h) const;

  /**
   * \fn Halfedge* get_mate( Halfedge* h ) const
   *
   * \brief  Returns the  mate  of  a given  half-edge  of this  PPS
   * underlying mesh.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the mate half-edge of a given half-edge.
   */
  Halfedge *get_mate(Halfedge *h) const;

  /**
   * \fn Halfedge* get_halfedge( Face* face ) const
   *
   * \brief Returns  the first half-edge of the  cycle of half-edges
   * of a given face of this PPS underlying mesh.
   *
   * \param face Pointer to a face of this PPS underlying mesh.
   *
   * \returns A pointer to the first half-edge of a given face.
   */
  Halfedge *get_halfedge(Face *face) const;

  /**
   * \fn Halfedge* get_halfedge( Vertex* vertex ) const
   *
   * \brief Returns one  half-edge with origin at a  given vertex of
   * this PPS underlying  mesh. It is assumed that  this method will
   * always return the same  half-edge (as many half-edges may share
   * the same origin vertex).
   *
   * \param vertex Pointer to a vertex of this PPS underlying mesh.
   *
   * \returns  A pointer  to one  half-edge with  origin at  a given
   * vertex.
   */
  Halfedge *get_halfedge(Vertex *vertex) const;

  /**
   * \fn unsigned get_degree( Vertex* vertex ) const
   *
   * \brief  Returns  the degree  of  a  given  vertex of  this  PPS
   * underlying mesh. The degree of  a vertex is the number of edges
   * incident to the vertex.
   *
   * \param vertex Pointer to a vertex of this PPS underlying mesh.
   *
   * \returns The degree of the given vertex.
   */
  unsigned get_degree(Vertex *vertex) const;

  /**
   * \fn Bezier* get_shape_function( Vertex* vertex ) const
   *
   * \brief  Returns  the shape  function  associated  with a  given
   * vertex of  this PPS  underlying mesh. The  shape function  is a
   * rectangular B&eacute;zier patch.
   *
   * \param vertex Pointer to a vertex of this PPS underlying mesh.
   *
   * \returns The shape function associated with the given vertex.
   */
  Bezier *get_shape_function(Vertex *vertex) const;

  /**
   * \fn void set_shape_function( Vertex* vertex, Bezier* patch )
   *
   * \brief Assigns a shape function  to a given vertex of this PPS
   * underlying mesh.
   *
   * \param vertex Pointer to a vertex of this PPS underlying mesh.
   * \param patch Pointer to a shape function.
   */
  void set_shape_function(Vertex *vertex, Bezier *patch);

  /**
   * \fn VertexIterator vertices_begin() const
   *
   * \brief Returns a vertex iterator set to the initial vertex of a
   * vertex sequence of this PPS underlying mesh.
   *
   * \return A vertex iterator set to the initial vertex of a vertex
   * sequence of this PPS underlying mesh.
   */
  VertexIterator vertices_begin() const;

  /**
   * \fn bool is_done( const VertexIterator& iterator ) const
   *
   * \brief Returns  a logic value  true if a given  vertex iterator
   * has reached the end of a vertex sequence of this PPS underlying
   * mesh; otherwise, it returns the logic value false.
   *
   * \param iterator A vertex iterator.
   *
   * \return  A logic  value true  if  a given  vertex iterator  has
   * reached the  end of  a vertex sequence  of this  PPS underlying
   * mesh; otherwise, it returns the logic value false.
   */
  bool is_done(const VertexIterator &iterator) const;

  /**
   * \fn void move_forward( VertexIterator& iterator ) const
   *
   * \brief Makes  the iterator point  to the vertex  succeeding its
   * current  vertex in  a vertex  sequence of  this  PPS underlying
   * mesh.
   *
   * \param iterator A reference to a vertex iterator.
   */
  void move_forward(VertexIterator &iterator) const;

  /**
   * \fn Vertex* get_vertex( const VertexIterator& iterator ) const
   *
   * \brief Returns  the current vertex  of a given  vertex iterator
   * for a vertex sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to a vertex iterator.
   *
   * \return  A pointer  to the  current  vertex of  a given  vertex
   * iterator for a vertex sequence of this PPS underlying mesh.
   */
  Vertex *get_vertex(const VertexIterator &iterator) const;

  /**
   * \fn EdgeIterator edges_begin() const
   *
   * \brief Returns an  edge iterator set to the  initial edge of an
   * edge sequence of this PPS underlying mesh.
   *
   * \return An  edge iterator  set to the  initial edge of  an edge
   * sequence of this PPS underlying mesh.
   */
  EdgeIterator edges_begin() const;

  /**
   * \fn bool is_done( const EdgeIterator& iterator ) const
   *
   * \brief Returns A logic value  true if a given edge iterator has
   * reached  the end  of an  edge sequence  of this  PPS underlying
   * mesh; otherwise, it returns the logic value false.
   *
   * \param iterator An edge iterator.
   *
   * \return A logic value true if a given edge iterator has reached
   * the  end of  an  edge  sequence of  this  PPS underlying  mesh;
   * otherwise, it returns the logic value false.
   */
  bool is_done(const EdgeIterator &iterator) const;

  /**
   * \fn void move_forward( EdgeIterator& iterator ) const
   *
   * \brief  Makes the  iterator point  to the  edge  succeeding its
   * current edge in an edge sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to an edge iterator.
   */
  void move_forward(EdgeIterator &iterator) const;

  /**
   * \fn Edge* get_edge( const EdgeIterator& iterator ) const
   *
   * \brief Returns the current edge of a given edge iterator for an
   * edge sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to an edge iterator.
   *
   * \return A pointer to the  current edge of a given edge iterator
   * for an edge sequence of this PPS underlying mesh.
   */
  Edge *get_edge(const EdgeIterator &iterator) const;

  /**
   * \fn FaceIterator faces_begin() const
   *
   * \brief Returns  a face  iterator set to  the initial face  of a
   * face sequence of this PPS underlying mesh.
   *
   * \return  A face  iterator set  to the  initial face  of  a face
   * sequence of this PPS underlying mesh.
   */
  FaceIterator faces_begin() const;

  /**
   * \fn bool is_done( const FaceIterator& iterator ) const
   *
   * \brief Returns A logic value  true if a given face iterator has
   * reached the end of a face sequence of this PPS underlying mesh;
   * otherwise, it returns the logic value false.
   *
   * \param iterator A face iterator.
   *
   * \return A logic value true if a given face iterator has reached
   * the  end  of a  face  sequence  of  this PPS  underlying  mesh;
   * otherwise, it returns the logic value false.
   */
  bool is_done(const FaceIterator &iterator) const;

  /**
   * \fn void move_forward( FaceIterator& iterator ) const
   *
   * \brief  Makes the  iterator point  to the  face  succeeding its
   * current face in a face sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to a face iterator.
   */
  void move_forward(FaceIterator &iterator) const;

  /**
   * \fn Face* get_face( const FaceIterator& iterator ) const
   *
   * \brief Returns the current face  of a given face iterator for a
   * face sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to a face iterator.
   *
   * \return A pointer to the  current face of a given face iterator
   * for a face sequence of this PPS underlying mesh.
   */
  Face *get_face(const FaceIterator &iterator) const;

private:
  /**
   * \fn void build_surface_evaluator()
   *
   * \brief Creates a Loop subdivision surface evaluator.
   */
  void build_surface_evaluator();

private:
  /**  Pointer to the Loop subdivision surface evaluator. */
  EvalStruct *_evstruct;
};

} // namespace ppsfromloop

/** @} */ // end of group class.
