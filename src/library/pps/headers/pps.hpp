#pragma once

/**
 * \file pps.hpp
 *
 * \brief  Definition of  the  abstract class  PPS  that represents  a
 * Parametric Pseudo Surface defined  on a generic underlying triangle
 * mesh.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 2.0
 * \date June 2009
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

#include <exception.hpp> // ERROR_UNLESS

#include <cmath> // std::abs, std::sqrt, std::acos, std::asin, std::cos, std::sin
#include <vector> // std::vector

/**
 * \defgroup PPSNameSpace Namespace pps.
 * @{
 */

/**
 * \namespace pps
 *
 * \brief The namespace pps  contains the definition and implementation
 * of all classes of  the open source library Parametric Pseudo-Surface
 * (PPS).
 *
 */

namespace pps {

/**
 * \class PPS
 *
 * \brief This class represents a Parametric Pseudo-Surface.
 *
 * The Mesh type is expected to  be a triangle mesh.  To use the PPS
 * class one  must implement  twenty-five pure virtual  methods. All
 * but one of  these methods ( "eval_surface" )  are related to mesh
 * topological  operations.  The  method  "eval_surface" computes  a
 * point on a generic parametric  patch associated with a given mesh
 * face.
 *
 * If one does not intend to use the PPS class with her/his own mesh
 * class, one can always use the default mesh class that accompanies
 * the PPS  class. For an example of  the use of the  PPS class with
 * the default mesh class, take a look at the classes PPSfromPNT and
 * PPSfromLOOP, both  of which accompanies  the PPS class  code. The
 * former  class  implements  the  "eval_surface"  method  as  a  PN
 * triangle patch, while the latter  implements the method as a Loop
 * surface patch.
 *
 */
template <typename Mesh> class PPS {
public:
  // ---------------------------------------------------------------
  //
  // Type name definitions for the generic mesh components.
  //
  // ---------------------------------------------------------------

  /**
   * \typedef Vertex
   *
   * \brief Definition of a type name for the mesh vertices.
   */
  using Vertex = typename Mesh::Vertex;

  /**
   * \using Halfedge
   *
   * \brief Definition of a type name for the mesh half-edges.
   */
  using Halfedge = typename Mesh::Halfedge;

  /**
   * \using Edge
   *
   * \brief Definition of a type name for the mesh edges.
   */
  using Edge = typename Mesh::Edge;

  /**
   * \using Face
   *
   * \brief Definition of a type name for the mesh faces.
   */
  using Face = typename Mesh::Face;

  /**
   * \using VertexIterator
   *
   * \brief Definition of a type name for the vertex iterators.
   */
  using VertexIterator = typename Mesh::VertexIterator;

  /**
   * \using EdgeIterator
   *
   * \brief Definition of a type name for the edge iterators.
   */
  using EdgeIterator = typename Mesh::EdgeIterator;

  /**
   * \using FaceIterator.
   *
   * \brief Definition of a type name for the face iterators.
   */
  using FaceIterator = typename Mesh::FaceIterator;

  // ---------------------------------------------------------------
  //
  // Public methods
  //
  // ---------------------------------------------------------------

  /**
   * \fn PPS( Mesh* mesh )
   *
   * \brief Creates an instance of this PPS class.
   *
   * \param mesh The address of an object of parameter type Mesh.
   *
   */
  PPS(Mesh *mesh) : _mesh(mesh), _MYPI(std::acos(-1)) {}

  /**
   * \fn virtual ~PPS()
   *
   * \brief Virtual destructor.
   */
  virtual ~PPS() {}

  /**
   * \fn void build()
   *
   * \brief Computes the parametrizations of this PPS, each of which
   * is associated  with one vertex of the  underlying triangle mesh
   * of the PPS.
   *
   */
  void build();

  /**
   * \fn void eval_pps( Face* face , double u , double v , double w , double& x
   * , double& y , double& z ) const
   *
   * \brief Computes a point on the image of this PPS. The resulting
   * point is the image of a point in a p-domain that corresponds to
   * a point  (through an  implicit homomorphism) in  a face  of the
   * underlying mesh.
   *
   * \param face Pointer to one face of the underlying mesh.
   * \param u First barycentric coordinate of a point on the face.
   * \param v Second barycentric coordinate of a point on the face.
   * \param w Third barycentric coordinate of a point on the face.
   * \param  x First  Cartesian coordinate  of  a point  on the  PPS
   * image.
   * \param  y Second  Cartesian coordinate  of a  point on  the PPS
   * image.
   * \param  z Third  Cartesian coordinate  of  a point  on the  PPS
   * image.
   *
   */
  void eval_pps(Face *face, double u, double v, double w, double &x, double &y,
                double &z) const;

  /**
   * \fn virtual void eval_surface( Face* face , double u , double v , double w
   * , double& x , double& y , double& z ) const = 0
   *
   * \brief Computes a point on the parametric patch associated with
   * a given face of the  PPS underlying triangle mesh.  This method
   * must be  implemented by the user  of the PPS  library (i.e., in
   * the  concrete PPS class  that inherits  from this  abstract PPS
   * class.
   *
   * \param face Pointer to one face of the underlying mesh.
   * \param u First barycentric coordinate of a point on the face.
   * \param v Second barycentric coordinate of a point on the face.
   * \param w Third barycentric coordinate of a point on the face.
   * \param x First Cartesian coordinate of a point on the PPS image.
   * \param y Second Cartesian coordinate of a point on the PPS image.
   * \param z Third Cartesian coordinate of a point on the PPS image.
   *
   */
  virtual void eval_surface(Face *face, double u, double v, double w, double &x,
                            double &y, double &z) const = 0;

  /**
   * \fn inline Mesh* get_mesh() const
   *
   * \brief Returns a pointer to the underlying mesh of this PPS.
   *
   * \return A pointer to the underlying mesh of this PPS.
   */
  inline Mesh *get_mesh() const { return _mesh; }

protected:
  // ---------------------------------------------------------------
  //
  // Protected methods
  //
  // ---------------------------------------------------------------

  /**
   * \fn void build_shape_functions()
   *
   * \brief  Creates   the  shape  functions   associated  with  the
   * p-domains.
   *
   */
  void build_shape_functions();

  /**
   * \fn void build_one_shape_function( Vertex* vertex )
   *
   * \brief  Creates  the shape  function  associated  with a  given
   * vertex of  this PPS underlying  mesh.  This method  samples the
   * p-domain associated with the  vertex and computes the points in
   * 3D that  are images  of the p-domain  points under  the generic
   * patch  defined  by   the  method  \link  pps::PPS::eval_surface
   * \endlink.  Those p-domain  points and their corresponding image
   * points are used to set up a linear system (the normal equations
   * of a least squares  problem), whose solution yields the control
   * points  of  the  shape  function.  This method  relies  on  the
   * implementation of \link pps::PPS::eval_surface \endlink.
   *
   * \param  vertex:  pointer  to  the vertex  associated  with  the
   * p-domain.
   *
   * \sa eval_surface()
   */
  void build_one_shape_function(Vertex *vertex);

  /**
   * \fn inline unsigned  get_shape_function_degree(  unsigned D  ) const
   *
   * \brief Returns  the bi-degree  of the rectangular  Bezier patch
   * defining  a shape  function.  This  bi-digree  is heuristically
   * defined as  the maximum between  the number seven and  one plus
   * the degree of the vertex associated with the shape function.
   *
   * \param D  The degree  of the vertex  associated with  the shape
   * function.
   *
   * \return The  bi-degree of the shape function  associated with a
   * vertex of a given degree.
   *
   */
  inline unsigned get_shape_function_degree(unsigned D) const {
    return (D > 6) ? 7 : (D + 1);
  }

  /**
   * \fn inline  unsigned get_number_of_parameter_points( unsigned D ) const
   *
   * \brief Returns the number of  points used to sample a p-domain.
   * This number is heuristically  defined as twice the bi-degree of
   * the shape function associated with the p-domain.
   *
   * \param D bi-degree  ( D , D ) of  the shape function associated
   * with a p-domain.
   *
   * \return The number of points used to sample a p-domain.
   *
   */
  inline unsigned get_number_of_parameter_points(unsigned D) const {
    return (D << 1);
  }

  /**
   * \fn inline unsigned get_degree( Halfedge* h ) const
   *
   * \brief  Returns the  degree of  the  origin vertex  of a  given
   * half-edge.
   *
   * \param h A pointer to a halfedge of this PPS underlying mesh.
   *
   * \returns The degree of the origin vertex of the given halfedge.
   *
   */
  inline unsigned get_degree(Halfedge *h) const {
    return get_degree(get_org(h));
  }

  /**
   * \fn bool find_triangle( Halfedge*& h , double x , double y , double& u ,
   * double& v , double& w ) const
   *
   * \brief Finds  the triangle of the canonical  triangulation of a
   * p-domain  that  contains a  given  point  (if  such a  triangle
   * exists). The p-domain is specified by a pointer to a half-edge,
   * i.   e., the  p-domain is  the one  associated with  the origin
   * vertex of the given half-edge.
   *
   * \param h  Pointer to  a half-edge. The  canonical triangulation
   * associated with the p-domain defined  by the origin vertex of h
   * is expected to contain the  given point. If so, the parameter h
   * will contain  the half-edge  whose origin vertex  is associated
   * with  the p-domain  corresponding to  the first  vertex  of the
   * triangle found by the method.
   * \param x First Cartesian coordinate of the given point.
   * \param y Second Cartesian coordinate of the given point.
   * \param u  First barycentric coordinate of the  given point with
   * respect to the triangle that contains it (if any).
   * \param v Second barycentric  coordinate of the given point with
   * respect to the triangle that contains it (if any).
   * \param w Third barycentric coordinate  of the given  point with
   * respect to the triangle that contains it (if any).
   *
   * \return True if a triangle is found and false otherwise.
   *
   */
  bool find_triangle(Halfedge *&h, double x, double y, double &u, double &v,
                     double &w) const;

  /**
   * \fn  void get_barycentric_coordinates(double x0 , double y0 , double x1 ,
   * double y1 , double x2 , double y2 , double xp , double yp , double& u ,
   * double& v , double& w ) const
   *
   * \brief Computes  the barycentric  coordinates of a  given point
   * (in Cartesian  coordinates) with  respect to a  given reference
   * triangle.
   *
   * \param x0 First Cartesian coordinate of the first vertex of the
   * reference triangle.
   * \param y0  Second Cartesian coordinate  of the first  vertex of
   * the reference triangle.
   * \param x1  First Cartesian coordinate  of the second  vertex of
   * the reference triangle.
   * \param y1  Second Cartesian coordinate of the  second vertex of
   * the reference triangle.
   * \param x2 First Cartesian coordinate of the third vertex of the
   * reference triangle.
   * \param y2  Second Cartesian coordinate  of the third  vertex of
   * the reference triangle.
   * \param xp First Cartesian coordinate of the point.
   * \param yp Second Cartesian coordinate of the point.
   * \param u First barycentric coordinate of the point.
   * \param v Second barycentric coordinate of the point.
   * \param w Third barycentric coordinate of the point.
   *
   */
  void get_barycentric_coordinates(double x0, double y0, double x1, double y1,
                                   double x2, double y2, double xp, double yp,
                                   double &u, double &v, double &w) const;

  /**
   * \fn void eval_pps( Halfedge* h , double u , double v , double w , double& x
   * , double& y , double& z ) const
   *
   * \brief Computes a point on the image of this PPS. The resulting
   * point is  the image of a  point in a p-domain.  The p-domain is
   * specified  by a half-edge,  whose origin  vertex is  the vertex
   * associated with  the p-domain. The  point is supposed to  be in
   * the face  that contains the  given half-edge in  the underlying
   * mesh of this  PPS. We are given the  barycentric coordinates of
   * the point, which  are the same coordinates of  the point in the
   * upper triangle of the canonical domain.
   *
   * \param  h  Pointer  to a  halfedge  of  the  face of  this  PPS
   * underlying mesh.
   * \param u First barycentric coordinate of point in the canonical
   * domain.
   * \param  v  Second barycentric  coordinate  of  a  point in  the
   * canonical p-domain.
   * \param  w  Third  barycentric  coordinate  of a  point  in  the
   * canonical p-domain.
   * \param  x First  Cartesian coordinate  of  a point  on the  PPS
   * image.
   * \param  y Second  Cartesian coordinate  of a  point on  the PPS
   * image.
   * \param z Third Cartesian coordinate of a point on the PPS image.
   */
  void eval_pps(Halfedge *h, double u, double v, double w, double &x, double &y,
                double &z) const;

  /**
   * \fn void select_pdomain( Face* face , double& u , double& v , double& w ,
   * Halfedge*&  h ) const
   *
   * \brief  Finds one  p-domain that  contains the  parameter point
   * that  is the  image of  a given  point in  a face  of  this PPS
   * underlying mesh.
   *
   * \param face Pointer to one face of the underlying mesh.
   * \param u First barycentric coordinate of a point on the face.
   * \param v Second barycentric coordinate of a point on the face.
   * \param w Third barycentric coordinate of a point on the face.
   * \param h A reference to a pointer to a half-edge.
   */
  void select_pdomain(Face *face, double &u, double &v, double &w,
                      Halfedge *&h) const;

  /**
   * \fn void from_barycentric_to_Cartesian( double u , double v , double w ,
   * double x0 , double y0 , double x1 , double y1 , double x2 , double y2 ,
   * double& x , double& y ) const
   *
   * \brief Converts  the barycentric  coordinates of a  point, with
   * respect to a reference triangle, to Cartesian coordinates.
   *
   * \param u First barycentric coordinate of the point.
   * \param v Second barycentric coordinate of the point.
   * \param w Third barycentric coordinate of the point.
   * \param x0 First Cartesian coordinate of the first vertex of the
   * reference triangle.
   * \param y0  Second Cartesian coordinate  of the first  vertex of
   * the reference triangle.
   * \param x1  First Cartesian coordinate  of the second  vertex of
   * the reference triangle.
   * \param y1  Second Cartesian coordinate of the  second vertex of
   * the reference triangle.
   * \param x2 First coordinate of the third vertex of the reference
   * triangle.
   * \param  y2  Second  coordinate  of  the  third  vertex  of  the
   * reference triangle.
   * \param x First Cartesian coordinate of the resulting point.
   * \param y Second Cartesian coordinate of the resulting point.
   *
   */
  void from_barycentric_to_Cartesian(double u, double v, double w, double x0,
                                     double y0, double x1, double y1, double x2,
                                     double y2, double &x, double &y) const;

  /**
   * \fn void compute_pdomain_contribution( Halfedge* h , double u , double v ,
   * double& weight , double& x , double& y , double& z ) const
   *
   * \brief  Computes   the  contribution  of   the  shape  function
   * associated  with  a  given  p-domain  at a  given  point.   The
   * p-domain is specified by  a given half-edge whose origin vertex
   * is the vertex associated with  the p-domain. The given point is
   * expected  to   be  in  the   canonical  triangulation  triangle
   * corresponding to the face that contains the given half-edge.
   *
   * \param h Pointer to a given half-edge.
   * \param u  First Cartesian  coordinate of the  given point  in a
   * p-domain.
   * \param v  Second Cartesian coordinate  of the given point  in a
   * p-domain.
   * \param weight accumulator of the sum of weights.
   * \param x  First Cartesian coordinate of  the point representing
   * the contribution of the shape function.
   * \param y Second Cartesian  coordinate of the point representing
   * the contribution of the shape function.
   * \param z  Third Cartesian coordinate of  the point representing
   * the contribution of the shape function.
   *
   */
  void compute_pdomain_contribution(Halfedge *h, double u, double v,
                                    double &weight, double &x, double &y,
                                    double &z) const;

  /**
   * \fn double gfunction( Halfedge* h , double x , double y , double& u ,
   * double& v ) const
   *
   * \brief  Computes the Cartesian  coordinates of  the image  of a
   * point  under  the transition  map  from  one  gluing domain  to
   * another.  The  source and target  p-domains are specified  by a
   * given half-edge. The source gluing domain is the one associated
   * with  the origin  vertex  of the  half-edge,  while the  target
   * gluing domain is the one associated with the destination vertex
   * of the given half-edge.
   *
   * \param h Pointer to a half-edge.
   * \param x First Cartesian coordinate of the given point.
   * \param y Second Cartesian coordinate of the given point.
   * \param u First  Cartesian coordinate of the image  of the given
   * point under the transition map.
   * \param v Second Cartesian coordinate of the image  of the given
   * point under the transition map.
   *
   */
  void gfunction(Halfedge *h, double x, double y, double &u, double &v) const;

  /**
   * \fn void to_canonical_domain( Halfedge* h , double x , double y , double& u
   * , double& v ) const
   *
   * \brief Maps a  point in a p-domain to a  point in the canonical
   * lens.  The p-domain is  specified by  a half-edge  whose origin
   * vertex is the one associated with the p-domain.
   *
   * \param h Pointer to a half-edge.
   * \param x First Cartesian coordinate of the given point.
   * \param y Second Cartesian coordinate of the given point.
   * \param u First  Cartesian coordinate of the image  of the given
   * point under the map.
   * \param v Second Cartesian coordinate of the image  of the given
   * point under the map.
   *
   */
  void to_canonical_domain(Halfedge *h, double x, double y, double &u,
                           double &v) const;

  /**
   * \fn void from_canonical_domain( Halfedge* h , double x , double y , double&
   * u , double& v ) const
   *
   * \brief Maps a  point in the canonical lens  to a p-domain.  The
   * p-domain is specified by a half-edge whose origin vertex is the
   * one associated with the p-domain.
   *
   * \param h Pointer to a half-edge.
   * \param x First Cartesian coordinate of the given point.
   * \param y Second Cartesian coordinate of the given point.
   * \param u First  Cartesian coordinate of the image  of the given
   * point under the map.
   * \param v Second Cartesian coordinate of the image  of the given
   * point under the map.
   *
   */
  void from_canonical_domain(Halfedge *h, double x, double y, double &u,
                             double &v) const;

  /**
   * \fn void rot_2d( double x , double y , double ang , double& u , double& v )
   * const
   *
   * \brief Computes the location of  a given point after a rotation
   * around the origin by a given angle.
   *
   * \param x First Cartesian coordinate of the given point.
   * \param y Second Cartesian coordinate of the given point.
   * \param ang A given rotation angle.
   * \param u First Cartesian coordinate of the point resulting from
   * the rotation of (x,y) by ang around the origin.
   * \param  v Second  Cartesian coordinate  of the  point resulting
   * from the rotation of (x,y) by ang around the origin.
   *
   */
  void rot_2d(double x, double y, double ang, double &u, double &v) const;

  /**
   * \fn  unsigned int get_id( Halfedge* h ) const
   *
   * \brief  Returns  the  identifier  of a  given  half-edge.   The
   * identifier of a half-edge is a  number ranging from 0 to N - 1,
   * where  N   is  the   degree  of  the   origin  vertex   of  the
   * half-edge. Each  such a number  represents the position  of one
   * half-edge  in  the   sequence  defined  by  a  counterclockwise
   * traversal of all half-edges with origin at the same vertex. The
   * first half-edge  visited in this traversal  gets the identifier
   * 0, the  second half-edge  gets the identifier  1, and so  on so
   * forth.
   *
   * \param h Pointer to a half-edge.
   *
   * \returns The identifier of the given half-edge.
   *
   */
  virtual unsigned int get_id(Halfedge *h) const;

  /**
   * \fn double weight_function( double x , double y , double R ) const
   *
   * \brief Evaluates the weight function at a given point.
   *
   * \param x First Cartesian coordinate of the given point.
   * \param y Second Cartesian coordinate of the given point.
   * \param R Radius of the circular support of the function.
   *
   * \return The function value at the given point.
   */
  double weight_function(double x, double y, double R) const;

  /**
   * \fn double eta_function( double s , double d1 , double d2 ) const
   *
   * \brief Evaluates the eta blending function at a given parameter
   * value.
   *
   * \param s Parameter value.
   * \param d1 Lower bound for the extent of the blending region.
   * \param d2 Upper bound for the extent of the blending region
   *
   * \return The function value at the given parameter value.
   */
  double eta_function(double s, double d1, double d2) const;

  // ---------------------------------------------------------------
  //
  // Abstract methods - must be implemented by the PPS library user.
  //
  // ---------------------------------------------------------------

  /**
   * \fn virtual bool mesh_has_boundary() const
   *
   * \brief  Determines if  the underlying  mesh of  this PPS  has a
   * non-empty boundary.
   *
   * \returns  The logic  value true  if  the mesh  had a  non-empty
   * boundary and the logic value false otherwise.
   *
   */
  virtual bool mesh_has_boundary() const = 0;

  /**
   * \fn virtual bool mesh_is_simplicial() const
   *
   * \brief  Determines if  the underlying  mesh  of this  PPS is  a
   * simplicial complex.
   *
   * \returns  The logic  value true  if  the mesh  is a  simplicial
   * complex and the logic value false otherwise.
   *
   */
  virtual bool mesh_is_simplicial() const = 0;

  /**
   * \fn virtual Vertex* get_org( Halfedge* h ) const
   *
   * \brief Returns the  origin vertex of a given  half-edge of this
   * PPS underlying mesh.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the origin vertex of a given half-edge.
   *
   */
  virtual Vertex *get_org(Halfedge *h) const = 0;

  /**
   * \fn virtual Vertex* get_dst( Halfedge* h ) const
   *
   * \brief Returns  the destination vertex of a  given half-edge of
   * this PPS underlying mesh.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns  A  pointer  to  the  destination vertex  of  a  given
   * half-edge.
   *
   */
  virtual Vertex *get_dst(Halfedge *h) const = 0;

  /**
   * \fn virtual Edge* get_edge( Halfedge* h ) const
   *
   * \brief  Returns  the  edge   a  given  half-edge  of  this  PPS
   * underlying mesh belongs to.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the edge a given half-edge belongs to.
   *
   */
  virtual Edge *get_edge(Halfedge *h) const = 0;

  /**
   * \fn virtual Face* get_face( Halfedge* h ) const
   *
   * \brief  Returns  the  face   a  given  half-edge  of  this  PPS
   * underlying mesh belongs to.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the face a given half-edge belongs to.
   *
   */
  virtual Face *get_face(Halfedge *h) const = 0;

  /**
   * \fn virtual Halfedge* get_prev( Halfedge* h ) const
   *
   * \brief Returns the half-edge that precedes a given half-edge of
   * this PPS underlying  mesh in the face cycle  of half-edges that
   * contains both half-edges.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns  A pointer  to  the half-edge  that  precedes a  given
   * half-edge in their common face half-edge cycle.
   *
   */
  virtual Halfedge *get_prev(Halfedge *h) const = 0;

  /**
   * \fn virtual Halfedge* get_next( Halfedge* h ) const
   *
   * \brief Returns the half-edge that succeeds a given half-edge of
   * this PPS underlying  mesh in the face cycle  of half-edges that
   * contains both half-edges.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns  A pointer  to  the half-edge  that  succeeds a  given
   * half-edge in their common face half-edge cycle.
   *
   */
  virtual Halfedge *get_next(Halfedge *h) const = 0;

  /**
   * \fn virtual Halfedge* get_mate( Halfedge* h ) const
   *
   * \brief  Returns the  mate  of  a given  half-edge  of this  PPS
   * underlying mesh.
   *
   * \param h Pointer to a half-edge of this PPS underlying mesh.
   *
   * \returns A pointer to the mate half-edge of a given half-edge.
   *
   */
  virtual Halfedge *get_mate(Halfedge *h) const = 0;

  /**
   * \fn virtual Halfedge* get_halfedge( Face* face ) const
   *
   * \brief Returns  the first half-edge of the  cycle of half-edges
   * of a given face of this PPS underlying mesh.
   *
   * \param face Pointer to a face of this PPS underlying mesh.
   *
   * \returns A pointer to the first half-edge of a given face.
   *
   */
  virtual Halfedge *get_halfedge(Face *face) const = 0;

  /**
   * \fn virtual Halfedge* get_halfedge( Vertex* vertex ) const
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
   *
   */
  virtual Halfedge *get_halfedge(Vertex *vertex) const = 0;

  /**
   * \fn virtual unsigned int get_degree( Vertex* vertex ) const
   *
   * \brief  Returns  the degree  of  a  given  vertex of  this  PPS
   * underlying mesh. The degree of  a vertex is the number of edges
   * incident to the vertex.
   *
   * \param vertex Pointer to a vertex of this PPS underlying mesh.
   *
   * \returns The degree of the given vertex.
   *
   */
  virtual unsigned int get_degree(Vertex *vertex) const = 0;

  /**
   * \fn virtual Bezier* get_shape_function( Vertex* vertex ) const
   *
   * \brief  Returns  the shape  function  associated  with a  given
   * vertex of  this PPS  underlying mesh. The  shape function  is a
   * rectangular B&eacute;zier patch.
   *
   * \param vertex Pointer to a vertex of this PPS underlying mesh.
   *
   * \returns The shape function associated with the given vertex.
   *
   */
  virtual Bezier *get_shape_function(Vertex *vertex) const = 0;

  /**
   * \fn virtual void set_shape_function( Vertex* vertex, Bezier* patch )
   *
   * \brief Assigns a  shape function to a given  vertex of this PPS
   * underlying mesh.
   *
   * \param vertex Pointer to a  vertex of this PPS underlying mesh.
   * \param patch Pointer to a shape function.
   *
   */
  virtual void set_shape_function(Vertex *vertex, Bezier *patch) = 0;

public:
  // ---------------------------------------------------------------
  //
  // Public abstract methods - must be implemented by the user.
  //
  // ---------------------------------------------------------------

  /**
   * \fn virtual VertexIterator vertices_begin() const
   *
   * \brief Returns a vertex iterator set to the initial vertex of a
   * vertex sequence of this PPS underlying mesh.
   *
   * \return A vertex iterator set to the initial vertex of a vertex
   * sequence of this PPS underlying mesh.
   *
   */
  virtual VertexIterator vertices_begin() const = 0;

  /**
   * \fn virtual bool is_done( const VertexIterator& iterator ) const
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
   *
   */
  virtual bool is_done(const VertexIterator &iterator) const = 0;

  /**
   * \fn virtual void move_forward( VertexIterator& iterator ) const
   *
   * \brief Makes  the iterator point  to the vertex  succeeding its
   * current  vertex in  a vertex  sequence of  this  PPS underlying
   * mesh.
   *
   * \param iterator A reference to a vertex iterator.
   *
   */
  virtual void move_forward(VertexIterator &iterator) const = 0;

  /**
   * \fn virtual Vertex* get_vertex( const VertexIterator& iterator ) const
   *
   * \brief Returns  the current vertex  of a given  vertex iterator
   * for a vertex sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to a vertex iterator.
   *
   * \return  A pointer  to the  current  vertex of  a given  vertex
   * iterator for a vertex sequence of this PPS underlying mesh.
   *
   */
  virtual Vertex *get_vertex(const VertexIterator &iterator) const = 0;

  /**
   * \fn virtual EdgeIterator edges_begin() const
   *
   * \brief Returns an  edge iterator set to the  initial edge of an
   * edge sequence of this PPS underlying mesh.
   *
   * \return An  edge iterator  set to the  initial edge of  an edge
   * sequence of this PPS underlying mesh.
   *
   */
  virtual EdgeIterator edges_begin() const = 0;

  /**
   * \fn virtual bool is_done( const EdgeIterator& iterator ) const
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
   *
   */
  virtual bool is_done(const EdgeIterator &iterator) const = 0;

  /**
   * \fn virtual void move_forward( EdgeIterator& iterator ) const
   *
   * \brief  Makes the  iterator point  to the  edge  succeeding its
   * current edge in an edge sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to an edge iterator.
   *
   */
  virtual void move_forward(EdgeIterator &iterator) const = 0;

  /**
   * \fn virtual Edge* get_edge( const EdgeIterator& iterator ) const
   *
   * \brief Returns the current edge of a given edge iterator for an
   * edge sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to an edge iterator.
   *
   * \return A pointer to the  current edge of a given edge iterator
   * for an edge sequence of this PPS underlying mesh.
   *
   */
  virtual Edge *get_edge(const EdgeIterator &iterator) const = 0;

  /**
   * \fn virtual FaceIterator faces_begin() const
   *
   * \brief Returns  a face  iterator set to  the initial face  of a
   * face sequence of this PPS underlying mesh.
   *
   * \return  A face  iterator set  to the  initial face  of  a face
   * sequence of this PPS underlying mesh.
   *
   */
  virtual FaceIterator faces_begin() const = 0;

  /**
   * \fn virtual bool is_done( const FaceIterator& iterator ) const
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
   *
   */
  virtual bool is_done(const FaceIterator &iterator) const = 0;

  /**
   * \fn virtual void move_forward( FaceIterator& iterator ) const
   *
   * \brief  Makes the  iterator point  to the  face  succeeding its
   * current face in a face sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to a face iterator.
   *
   */
  virtual void move_forward(FaceIterator &iterator) const = 0;

  /**
   * \fn virtual Face* get_face( const FaceIterator& iterator ) const
   *
   * \brief Returns the current face  of a given face iterator for a
   * face sequence of this PPS underlying mesh.
   *
   * \param iterator A reference to a face iterator.
   *
   * \return A pointer to the  current face of a given face iterator
   * for a face sequence of this PPS underlying mesh.
   *
   */
  virtual Face *get_face(const FaceIterator &iterator) const = 0;

protected:
  // ---------------------------------------------------------------
  //
  // Protected data members
  //
  // ---------------------------------------------------------------

  Mesh *_mesh; ///< A pointer to a triangle mesh of type Mesh.

  const double _MYPI; ///< The constant PI.
};

// -----------------------------------------------------------------
//
// Implementation of non-virtual public methods
//
// -----------------------------------------------------------------

/**
 * \fn void PPS< Mesh >::eval_pps( Face* face , double u , double v , double w ,
 * double& x , double& y , double& z ) const
 *
 * \brief Computes a  point on the image of  this PPS. The resulting
 * point is the image of a point in a p-domain that corresponds to a
 * point  (through  an  implicit  homomorphism)  in a  face  of  the
 * underlying mesh.
 *
 * \param face Pointer to one face of the underlying mesh.
 * \param u First barycentric coordinate of a point on the face.
 * \param v Second barycentric coordinate of a point on the face.
 * \param w Third barycentric coordinate of a point on the face.
 * \param x First Cartesian coordinate of a point on the PPS image.
 * \param y Second Cartesian coordinate of a point on the PPS image.
 * \param z Third Cartesian coordinate of a point on the PPS image.
 *
 */
template <typename Mesh>
void PPS<Mesh>::eval_pps(Face *face, double u, double v, double w, double &x,
                         double &y, double &z) const {
  /*
   * The barycentric coordinates must define a point in the face.
   */
  ERROR_UNLESS((u >= 0) && (u <= 1), "Invalid barycentric coordinates");
  ERROR_UNLESS((v >= 0) && (v <= 1), "Invalid barycentric coordinates");
  ERROR_UNLESS((w >= 0) && (w <= 1), "Invalid barycentric coordinates");
  ERROR_UNLESS(std::abs(1 - (u + v + w)) <= 1e-15,
               "Invalid barycentric coordinates");

  /*
   * We choose the  p-domain whose distance from the  given point is
   * the smallest.
   */

  Halfedge *h;
  double uaux = u;
  double vaux = v;
  double waux = w;

  select_pdomain(face, uaux, vaux, waux, h);

  /*
   * Evaluate the PPS using the chosen p-domain.
   */
  eval_pps(h, uaux, vaux, waux, x, y, z);

  return;
}

// -----------------------------------------------------------------
//
// Implementation of non-virtual protected methods.
//
// -----------------------------------------------------------------

/**
 * \fn void PPS< Mesh >::build()
 *
 * \brief Computes  the parametrizations of this PPS,  each of which
 * is associated with one vertex  of the underlying triangle mesh of
 * the PPS.
 *
 */
template <typename Mesh> void PPS<Mesh>::build() {
  /*
   * The given underlying mesh must be a triangle mesh.
   */
  ERROR_UNLESS(mesh_is_simplicial(), "Expected a simplicial mesh");

  /*
   * The given mesh must have an empty boundary.
   */
  ERROR_UNLESS(!mesh_has_boundary(), "Expected a mesh with an empty boundary");

  /*
   * Compute the shape functions of this PPS.
   */
  build_shape_functions();
}

/**
 * \fn void PPS< Mesh >::build_shape_functions()
 *
 * \brief   Creates  the   shape  functions   associated   with  the
 * p-domains.  This method  relies  on the  implementation of  \link
 * eval_surface() \endlink.
 *
 * \sa eval_surface()
 *
 */
template <typename Mesh> void PPS<Mesh>::build_shape_functions() {
  /*
   * For  each  vertex  "v"  of  the  underlying  mesh,  sample  the
   * P-polygon   associated   with   it,   and  then   compute   the
   * corresponding  points in  the generic  parametric  patch. These
   * points are then used to compute the control points of the shape
   * functions.
   */
  for (VertexIterator vi = vertices_begin(); !is_done(vi); move_forward(vi)) {
    /*
     * Get the current vertex.
     */
    Vertex *vertex = get_vertex(vi);

    /*
     * Compute  the  shape  function  corresponding to  the  current
     * vertex.
     */
    build_one_shape_function(vertex);
  }

  return;
}

/**
 * \fn void PPS< Mesh >::build_one_shape_function( Vertex* vertex )
 *
 * \brief Creates the shape  function associated with a given vertex
 * of this  PPS underlying mesh.   This method samples  the p-domain
 * associated with the vertex and computes the points in 3D that are
 * images of the p-domain points  under the generic patch defined by
 * the  method  \link  PPS<  Mesh >::eval_surface  \endlink.   Those
 * p-domain points and their  corresponding image points are used to
 * set up a  linear system (the normal equations  of a least squares
 * problem), whose  solution yields the control points  of the shape
 * function.
 *
 * \param  vertex:  pointer  to   the  vertex  associated  with  the
 * p-domain.
 *
 */
template <typename Mesh>
void PPS<Mesh>::build_one_shape_function(Vertex *vertex) {
#ifdef DEBUGMODE
  /**
   * \pre{ Pointer vertex cannot be null. }
   */

  ERROR_UNLESS(vertex != 0, "Expected a non-null pointer");
#endif

  /*
   * Get one halfedge with origin at vertex \var vertex.
   */
  Halfedge *h = get_halfedge(vertex);

#ifdef DEBUGMODE
  /**
   * \pre{ Pointer h cannot be null. }
   */

  ERROR_UNLESS(h != 0, "Expected a non-null pointer");
#endif

  unsigned int nu = get_degree(h);

#ifdef DEBUGMODE
  /**
   * \pre{ The degree of a vertex must be at least 3. }
   */

  ERROR_UNLESS(nu >= 3, "Vertex degree must be at least 3");
#endif

  /*
   * Generate a  rectangular grid with (N  + 1) * (N  + 1) parameter
   * points.
   */

  const unsigned int D = get_shape_function_degree(nu);
  const unsigned int N = get_number_of_parameter_points(D);

  const double R = std::cos(_MYPI / nu);

  /*
   * Coordinates of the lower leftmost point of the grid.
   */
  const double x0 = -R;
  const double y0 = -R;

  /*
   * Spacing  between  two  consecutives  points  in both  X  and  Y
   * directions.
   */
  const double dd = (2 * R) / N;

  /*
   * Initialize the Y direction spacing counter.
   */
  double dy = 0;

  std::vector<double> param_pts; // Array of parameter points.
  std::vector<double> patch_pts; // Array of surface points.

  for (unsigned int j = 0; j <= N; j++) {
    /* Y coordinate of the point. */
    double y = y0 + dy;

    /*
     * Initialize the X direction spacing counter.
     */
    double dx = 0;

    for (unsigned int i = 0; i <= N; i++) {
      /* X coordinate of the point. */
      double x = x0 + dx;

      /*
       * If  this  parameter  point  is  inside  the  P-polygon  the
       * P-domain is  inscribed in, then  we can compute a  point on
       * the given generic patch.
       */

      double u;
      double v;
      double w;
      Halfedge *haux = h;
      if (find_triangle(haux, x, y, u, v, w)) {
        /*
         * Store the coordinates of the parameter point.
         */
        param_pts.push_back(x);
        param_pts.push_back(y);

        /*
         * Compute  the corresponding point  in the  generic surface
         * patch.
         */
#ifdef DEBUGMODE
        ERROR_UNLESS(haux != 0, "Expected a non-null pointer");
#endif

        /* Get the face the half-edge haux belongs to. */
        Face *face = get_face(haux);

#ifdef DEBUGMODE
        ERROR_UNLESS(face != 0, "Expected a non-null pointer");
#endif

        double pt[3];

        if (haux == get_halfedge(face)) {
          eval_surface(face, u, v, w, pt[0], pt[1], pt[2]);
        } else if (get_next(haux) == get_halfedge(face)) {
          eval_surface(face, v, w, u, pt[0], pt[1], pt[2]);
        } else {
          eval_surface(face, w, u, v, pt[0], pt[1], pt[2]);
        }

        /*
         * Store the coordinates of the point on the generic surface
         * patch.
         */
        patch_pts.push_back(pt[0]);
        patch_pts.push_back(pt[1]);
        patch_pts.push_back(pt[2]);
      }

      /* Increment X direction spacing counter. */
      dx += dd;
    }

    /* Increment Y direction spacing counter. */
    dy += dd;
  }

  /* Creates the shape function. */
  set_shape_function(vertex,
                     new Bezier(&param_pts[0], &patch_pts[0],
                                param_pts.size() >> 1, D, D, -R, -R, R, R));

  return;
}

/**
 * \fn bool PPS< Mesh >::find_triangle( Halfedge*& h , double x , double y ,
 * double& u , double& v , double& w ) const
 *
 * \brief  Finds the triangle  of the  canonical triangulation  of a
 * p-domain  that  contains  a  given  point  (if  such  a  triangle
 * exists). The p-domain  is specified by a pointer  to a half-edge,
 * i.  e., the p-domain is the one associated with the origin vertex
 * of the given half-edge.
 *
 * \param  h Pointer  to  a half-edge.  The canonical  triangulation
 * associated with the p-domain defined by the origin vertex of h is
 * expected to contain the given  point. If so, the parameter h will
 * contain the half-edge whose  origin vertex is associated with the
 * p-domain corresponding to the  first vertex of the triangle found
 * by the method.
 * \param x First Cartesian coordinate of the given point.
 * \param y Second Cartesian coordinate of the given point.
 * \param  u First barycentric  coordinate of  the given  point with
 * respect to the triangle that contains it (if any).
 * \param v  Second barycentric coordinate  of the given  point with
 * respect to the triangle that contains it (if any).
 * \param  w Third barycentric  coordinate of  the given  point with
 * respect to the triangle that contains it (if any).
 *
 * \return True if a triangle is found and false otherwise.
 *
 */
template <typename Mesh>
bool PPS<Mesh>::find_triangle(Halfedge *&h, double x, double y, double &u,
                              double &v, double &w) const {
  /*
   * Get the degree of the vertex associated with the P-polygon.
   */
  unsigned int nu = get_degree(h);

  /*
   * Loop over  all triangles of the canonical  triangulation of the
   * P-polygon associated with the  origin vertex of halfedge h. For
   * each triangle, checks if the  given point belongs to it. If so,
   * compute the  barycentric coordinates of the  point with respect
   * to the triangle.
   */
  Halfedge *haux = h;
  do {

#ifdef DEBUGMODE
    ERROR_UNLESS(haux != 0, "Expected a non-null pointer");
#endif

    /*
     * Get the identifier of the current half-edge.
     */
    unsigned int id = get_id(haux);

    /*
     * Compute the vertices of the P-polygon canonical triangulation
     * triangle associated with the face that contains the half-edge
     * haux.
     */
    const double ang = 2 * (_MYPI / nu);

    double x1 = std::cos(id * ang);
    double y1 = std::sin(id * ang);
    double x2 = std::cos((id + 1) * ang);
    double y2 = std::sin((id + 1) * ang);

    /*
     * Compute the  barycentric coordinates of the  given point with
     * respect to the triangle  given by the vertices ( 0 ,  0 ) , (
     * x1 , y1 ) , and ( x2 , y2 ).
     */
    get_barycentric_coordinates(0., 0., x1, y1, x2, y2, x, y, u, v, w);

    /*
     * If all  barycentric coordinates are equal to  or greater than
     * zero,  then the  triangle contains  the given  point  and the
     * search  ends. Otherwise,  keep  looking for  a triangle  that
     * contains the point.
     */
    if ((u >= 0) && (u <= 1) && (v >= 0) && (v <= 1) && (w >= 0) && (w <= 1)) {
      h = haux;
      return true;
    } else {
      haux = get_mate(get_prev(haux));
    }

  } while (haux != h);

  /*
   * If the code  reached this point, then no  triangle contains the
   * given point.
   */
  return false;
}

/**
 * \fn void PPS< Mesh >::get_barycentric_coordinates(double x0 , double y0 ,
 * double x1 , double y1 , double x2 , double y2 , double xp , double yp ,
 * double& u , double& v , double& w ) const
 *
 * \brief Computes the barycentric  coordinates of a given point (in
 * Cartesian  coordinates)   with  respect  to   a  given  reference
 * triangle.
 *
 * \param x0 First  Cartesian coordinate of the first  vertex of the
 * reference triangle.
 * \param y0 Second Cartesian coordinate  of the first vertex of the
 * reference triangle.
 * \param x1 First Cartesian coordinate  of the second vertex of the
 * reference triangle.
 * \param y1 Second Cartesian coordinate of the second vertex of the
 * reference triangle.
 * \param x2 First  Cartesian coordinate of the third  vertex of the
 * reference triangle.
 * \param y2 Second Cartesian coordinate  of the third vertex of the
 * reference triangle.
 * \param xp First Cartesian coordinate of the point.
 * \param yp Second Cartesian coordinate of the point.
 * \param u First barycentric coordinate of the point.
 * \param v Second barycentric coordinate of the point.
 * \param w Third barycentric coordinate of the point.
 *
 */
template <typename Mesh>
void PPS<Mesh>::get_barycentric_coordinates(double x0, double y0, double x1,
                                            double y1, double x2, double y2,
                                            double xp, double yp, double &u,
                                            double &v, double &w) const {
  /*
   * Compute the determinant.
   */
  double dd =
      (x1 * y0) - (x2 * y0) - (x0 * y1) + (x2 * y1) + (x0 * y2) - (x1 * y2);

  /*
   * The determinant cannot be zero.
   */
  ERROR_UNLESS(std::abs(dd) > 1e-16, "Expected a non-zero determinant");

  /*
   * Compute the barycentric coordinates.
   */
  u = (x2 * y1) - (xp * y1) - (x1 * y2) + (xp * y2) + (x1 * yp) - (x2 * yp);

  u /= dd;

  v = (xp * y0) - (x2 * y0) + (x0 * y2) - (xp * y2) - (x0 * yp) + (x2 * yp);

  v /= dd;

  if (std::abs(u) < 1e-14) {
    u = 0;
  } else if (std::abs(1 - u) < 1e-14) {
    u = 1;
  }

  if (std::abs(v) < 1e-14) {
    v = 0;
  } else if (std::abs(1 - v) < 1e-14) {
    v = 1;
  }

  w = 1 - u - v;

  if (std::abs(w) < 1e-14) {
    w = 0;
  } else if (std::abs(1 - w) < 1e-14) {
    w = 1;
  }

  return;
}

/**
 * \fn void PPS< Mesh >::eval_pps( Halfedge* h , double u , double v , double w
 * , double& x , double& y , double& z ) const
 *
 * \brief Computes a  point on the image of  this PPS. The resulting
 * point is  the image  of a  point in a  p-domain. The  p-domain is
 * specified  by a  half-edge,  whose origin  vertex  is the  vertex
 * associated with the p-domain. The  point is supposed to be in the
 * face that contains the given  half-edge in the underlying mesh of
 * this PPS. We are given  the barycentric coordinates of the point,
 * which are the same coordinates of the point in the upper triangle
 * of the canonical domain.
 *
 * \param h Pointer to a halfedge of the face of this PPS underlying
 * mesh.
 * \param u  First barycentric coordinate of point  in the canonical
 * domain.
 * \param  v  Second  barycentric  coordinate  of  a  point  in  the
 * canonical p-domain.
 * \param w Third barycentric coordinate of a point in the canonical
 * p-domain.
 * \param x First Cartesian coordinate of a point on the PPS image.
 * \param y Second Cartesian coordinate of a point on the PPS image.
 * \param z Third Cartesian coordinate of a point on the PPS image.
 *
 */
template <typename Mesh>
void PPS<Mesh>::eval_pps(Halfedge *h, double u, double v, double w, double &x,
                         double &y, double &z) const {
#ifdef DEBUGMODE
  /**
   * \pre{ Pointer h cannot be null. }
   */

  ERROR_UNLESS(h != 0, "Expected a non-null pointer");
#endif

  /*
   * Computes  the Cartesian  coordinates  of the  given point  with
   * respect  to   the  upper  triangle  of   the  canonical  domain
   * (quadrilateral).
   */
  double xc;
  double yc;
  from_barycentric_to_Cartesian(u, v, w, 0., 0., 1., 0., 0.5,
                                0.5 * std::sqrt(3.), xc, yc);

  /*
   * Map  the  point  from  the  canonical domain  to  the  p-domain
   * associated with the origin vertex of h.
   */
  double xr;
  double yr;
  from_canonical_domain(h, xc, yc, xr, yr);

  /* Get the degree of the origin vertex of the given half-edge. */
  unsigned int nu = get_degree(h);

  /*
   * Compute the  radius of the p-domain associated  with the origin
   * vertex of the given half-edge.
   */
  const double R = std::cos(_MYPI / nu);

  double ll = (xr * xr) + (yr * yr);

  /*
   * The given point must belong to the p-domain.
   */
  ERROR_UNLESS(ll < R * R, "Given point does not belong to p-domain");

  /*
   * Initialize the accumulator of the sum of weight function values
   * and  the  accumulator of  the  contribution  of p-domain  shape
   * functions.
   */
  double sw = 0;
  double sf[3] = {0., 0., 0.};

  /*
   * Computes the  contribution of the p-domain  associated with the
   * origin vertex of h.
   */

  /* Compute the value of the weight function. */
  sw = weight_function(xr, yr, R);

  /*
   * The weight of the given point can never be zero.
   */
  ERROR_UNLESS(sw >= 1e-16, "The weight of a point is expected to be positive");

  /** Compute the value of the shape function. */
  Vertex *vertex = get_org(h);

#ifdef DEBUGMODE
  ERROR_UNLESS(vertex != 0, "Expected a non-null pointer");
#endif

  /*
   * Get the shape function associated with the origin vertex of the
   * given half-edge.
   */
  Bezier *patch = get_shape_function(vertex);

#ifdef DEBUGMODE
  ERROR_UNLESS(patch != 0, "Expected a non-null pointer");
#endif

  /*
   * Compute  a point  on  the patch  corresponding  to the  shape
   * function.
   */
  patch->point(xr, yr, sf[0], sf[1], sf[2]);

  sf[0] = sw * sf[0];
  sf[1] = sw * sf[1];
  sf[2] = sw * sf[2];

  /*
   * Compute the contribution of the other two p-domains (if any).
   */
  if ((u != 1) && (v != 1) && (w != 1)) {
    /*
     * If the point  lies on an edge, then we  can consider only one
     * more p-domain.
     */
    if (w == 0) {
      double ur;
      double vr;
      gfunction(h, xr, yr, ur, vr);

      compute_pdomain_contribution(get_mate(h), ur, vr, sw, sf[0], sf[1],
                                   sf[2]);
    } else if (v == 0) {
      double ur;
      double vr;
      Halfedge *h2 = get_mate(get_prev(h));

      gfunction(h2, xr, yr, ur, vr);

      compute_pdomain_contribution(get_mate(h2), ur, vr, sw, sf[0], sf[1],
                                   sf[2]);
    } else {
      double ur;
      double vr;
      gfunction(h, xr, yr, ur, vr);

      compute_pdomain_contribution(get_mate(h), ur, vr, sw, sf[0], sf[1],
                                   sf[2]);

      Halfedge *h2 = get_mate(get_prev(h));

      gfunction(h2, xr, yr, ur, vr);

      compute_pdomain_contribution(get_mate(h2), ur, vr, sw, sf[0], sf[1],
                                   sf[2]);
    }
  }

  /*
   * Compute the  point on the PPS after  weighting the contribution
   * of the three p-domains.
   */
  x = sf[0] / sw;
  y = sf[1] / sw;
  z = sf[2] / sw;

  return;
}

/**
 * \fn void PPS< Mesh >::select_pdomain( Face* face , double& u , double& v ,
 * double& w , Halfedge*&  h ) const
 *
 * \brief Finds one p-domain  that contains the parameter point that
 * is the  image of a given point  in a face of  this PPS underlying
 * mesh.
 *
 * \param face Pointer to one face of the underlying mesh.
 * \param u First barycentric coordinate of a point on the face.
 * \param v Second barycentric coordinate of a point on the face.
 * \param w Third barycentric coordinate of a point on the face.
 * \param h A reference to a pointer to a half-edge.
 */
template <typename Mesh>
void PPS<Mesh>::select_pdomain(Face *face, double &u, double &v, double &w,
                               Halfedge *&h) const {
#ifdef DEBUGMODE
  /**
   * \pre{ Pointer face cannot be null. }
   */

  ERROR_UNLESS(face != 0, "Expected a non-null pointer");
#endif

  /*
   * Find out which P-domain  contains the given point after mapping
   * the point to the canonical quadrilateral.
   */
  double xc;
  double yc;
  from_barycentric_to_Cartesian(u, v, w, 0., 0., 1., 0., 0.5,
                                0.5 * std::sqrt(3.), xc, yc);

  /*
   * We choose the  p-domain whose distance from the  given point is
   * the smallest.
   */
  double l1 = (xc * xc) + (yc * yc);

  double x2 = xc - 1;
  double l2 = (x2 * x2) + (yc * yc);

  double x3 = xc - 0.5;
  double y3 = yc - (0.5 * std::sqrt(3.));

  double l3 = (x3 * x3) + (y3 * y3);

  if ((l1 <= l2) && (l1 <= l3)) {
    /*
     * We pick the p-domain associated with the origin vertex of the
     * first half-edge of the given face.
     */
    h = get_halfedge(face);
  } else if (l2 <= l3) {
    /*
     * We pick the p-domain associated with the origin vertex of the
     * second half-edge of the given face.
     */
    h = get_next(get_halfedge(face));

    double uaux = u;
    u = v;
    v = w;
    w = uaux;
  } else {
    /*
     * We pick the p-domain associated with the origin vertex of the
     * third half-edge of the given face.
     */
    h = get_prev(get_halfedge(face));

    double uaux = u;
    u = w;
    w = v;
    v = uaux;
  }

  return;
}

/**
 * \fn void PPS< Mesh >::from_barycentric_to_Cartesian( double u , double v ,
 * double w , double x0 , double y0 , double x1 , double y1 , double x2 , double
 * y2 , double& x , double& y ) const
 *
 * \brief  Converts the  barycentric  coordinates of  a point,  with
 * respect to a reference triangle, to Cartesian coordinates.
 *
 * \param u First barycentric coordinate of the point.
 * \param v Second barycentric coordinate of the point.
 * \param w Third barycentric coordinate of the point.
 * \param x0 First  Cartesian coordinate of the first  vertex of the
 * reference triangle.
 * \param y0 Second Cartesian coordinate  of the first vertex of the
 * reference triangle.
 * \param x1 First Cartesian coordinate  of the second vertex of the
 * reference triangle.
 * \param y1 Second Cartesian coordinate of the second vertex of the
 * reference triangle.
 * \param x2 First  coordinate of the third vertex  of the reference
 * triangle.
 * \param y2 Second coordinate of  the third vertex of the reference
 * triangle.
 * \param x First Cartesian coordinate of the resulting point.
 * \param y Second Cartesian coordinate of the resulting point.
 *
 */
template <typename Mesh>
void PPS<Mesh>::from_barycentric_to_Cartesian(double u, double v, double w,
                                              double x0, double y0, double x1,
                                              double y1, double x2, double y2,
                                              double &x, double &y) const {
  x = (u * x0) + (v * x1) + (w * x2);
  y = (u * y0) + (v * y1) + (w * y2);
}

/**
 * \fn void PPS< Mesh >::compute_pdomain_contribution( Halfedge* h , double u ,
 * double v , double& weight , double& x , double& y , double& z ) const
 *
 * \brief Computes the contribution of the shape function associated
 * with  a  given  p-domain  at  a given  point.   The  p-domain  is
 * specified by a given half-edge  whose origin vertex is the vertex
 * associated with the  p-domain. The given point is  expected to be
 * in the canonical triangulation triangle corresponding to the face
 * that contains the given half-edge.
 *
 * \param h Pointer to a  given half-edge.
 * \param  u First  Cartesian coordinate  of  the given  point in  a
 * p-domain.
 * \param  v Second  Cartesian coordinate  of the  given point  in a
 * p-domain.
 * \param weight accumulator of the sum of weights.
 * \param x First Cartesian coordinate of the point representing the
 * contribution of the shape function.
 * \param y  Second Cartesian  coordinate of the  point representing
 * the contribution of the shape function.
 * \param z Third Cartesian coordinate of the point representing the
 * contribution of the shape function.
 *
 */
template <typename Mesh>
void PPS<Mesh>::compute_pdomain_contribution(Halfedge *h, double u, double v,
                                             double &weight, double &x,
                                             double &y, double &z) const {
  /*
   * Get the degree of the origin vertex of the given half-edge.
   */
  unsigned int nu = get_degree(h);

  /*
   * Get the  radius of the  boundary circumference of  the p-domain
   * associated with the origin vertex of the given half-edge.
   */
  const double R = std::cos(_MYPI / nu);

  /*
   * Compute the squared  distance of the point (u,v)  to the origin
   * of the local coordinate  system of the p-domain associated with
   * the given half-edge.
   */
  double ll = u * u + v * v;

  /*
   * If the  point (u,v)  is inside the  p-domain, then  compute the
   * contribution of  the shape function associated  with the origin
   * vertex of the given half-edge.
   */
  if (ll < R * R) {
    /*
     * Compute the weight of the contribution.
     */
    double w = weight_function(u, v, R);
    weight += w;

    /*
     * Get the  shape function associated with the  origin vertex of
     * the given half-edge.
     */
    Bezier *patch = get_shape_function(get_org(h));

#ifdef DEBUGMODE
    ERROR_UNLESS(patch != 0, "Expected a non-null pointer");
#endif

    /*
     * Compute  a point  on  the patch  corresponding  to the  shape
     * function.
     */
    double fx, fy, fz;
    patch->point(u, v, fx, fy, fz);

    /*
     * Weight the contribution of the shape function.
     */
    x += w * fx;
    y += w * fy;
    z += w * fz;
  }

  return;
}

/**
 * \fn double PPS< Mesh >::gfunction( Halfedge* h , double x , double y ,
 * double& u , double& v ) const
 *
 * \brief Computes the Cartesian coordinates of the image of a point
 * under the transition map from  one gluing domain to another.  The
 * source   and  target   p-domains   are  specified   by  a   given
 * half-edge. The  source gluing domain  is the one  associated with
 * the  origin vertex  of  the half-edge,  while  the target  gluing
 * domain is the  one associated with the destination  vertex of the
 * given half-edge.
 *
 * \param h Pointer to a half-edge.
 * \param x First Cartesian coordinate of the given point.
 * \param y Second Cartesian coordinate of the given point.
 * \param u  First Cartesian  coordinate of the  image of  the given
 * point under the transition map.
 * \param v  Second Cartesian coordinate  of the image of  the given
 * point under the transition map.
 *
 */
template <typename Mesh>
void PPS<Mesh>::gfunction(Halfedge *h, double x, double y, double &u,
                          double &v) const {

#ifdef DEBUGMODE
  ERROR_UNLESS(h != 0, "Expected a non-null pointer");
#endif

  /*
   * Map  the point  in the  source gluing  domain to  the canonical
   * lens.
   */
  double utemp;
  double vtemp;
  to_canonical_domain(h, x, y, utemp, vtemp);

  /*
   * Apply the reflection in the canonical quadrilateral.
   */
  utemp = 1 - utemp;
  vtemp = -vtemp;

  /*
   * Map the point in the canonical lens to the target p-domain.
   */
  from_canonical_domain(get_mate(h), utemp, vtemp, u, v);
}

/**
 * \fn void PPS< Mesh >::to_canonical_domain( Halfedge* h , double x , double y
 * , double& u , double& v ) const
 *
 * \brief Maps  a point in  a p-domain to  a point in  the canonical
 * domain.   The p-domain  is specified  by a  half-edge  whose origin
 * vertex is the one associated with the p-domain.
 *
 * \param h Pointer to a half-edge.
 * \param x First Cartesian coordinate of the given point.
 * \param y Second Cartesian coordinate of the given point.
 * \param u  First Cartesian  coordinate of the  image of  the given
 * point under the map.
 * \param v  Second Cartesian coordinate  of the image of  the given
 * point under the map.
 *
 */
template <typename Mesh>
void PPS<Mesh>::to_canonical_domain(Halfedge *h, double x, double y, double &u,
                                    double &v) const {
  /*
   * Get the degree of the origin vertex of the given half-edge.
   */
  const unsigned int nu = get_degree(h);

  /*
   * Get the rotation angle associated with the given half-edge.
   */
  const double au = (2 * _MYPI) / nu;

  /*
   * Get the identifier of the given half-edge.
   */
  const unsigned i = get_id(h);

  /*
   * Perform  a 2D  rotation by  -i *  au around  the origin  of the
   * p-domain  associated  with  the  origin  vertex  of  the  given
   * half-edge.
   */
  rot_2d(x, y, -(i * au), u, v);

  /*
   * Change from Cartesian to polar coordinates.
   */
  double rr = std::sqrt((u * u) + (v * v));

  if (std::abs(rr) > 1e-16) {
    double aa;

    if (u < 0) {
      if (v >= 0) {
        aa = _MYPI - std::acos(std::abs(u) / rr);
      } else {
        aa = (2 * _MYPI) - std::acos(std::abs(u) / rr);
      }
    } else {
      if (v >= 0) {
        aa = std::asin(v / rr);
      } else {
        aa = -std::asin(std::abs(v) / rr);
      }
    }

    aa *= (nu / 6.);

    rr *= (std::cos(_MYPI / 6.) / std::cos(_MYPI / nu));

    u = rr * std::cos(aa);
    v = rr * std::sin(aa);
  } else {
    u = v = 0.;
  }

  if (std::abs(u) < 1e-15)
    u = 0.;

  if (std::abs(v) < 1e-15)
    v = 0.;

  if (std::abs(1 - u) < 1e-15)
    u = 1.;

  if (std::abs(1 - v) < 1e-15)
    v = 1.;

  return;
}

/**
 * \fn void PPS< Mesh >::from_canonical_domain( Halfedge* h , double x , double
 * y , double& u , double& v ) const
 *
 * \brief Maps  a point  in the canonical  lens to a  p-domain.  The
 * p-domain is specified  by a half-edge whose origin  vertex is the
 * one associated with the p-domain.
 *
 * \param h Pointer to a half-edge.
 * \param x First Cartesian coordinate of the given point.
 * \param y Second Cartesian coordinate of the given point.
 * \param u  First Cartesian  coordinate of the  image of  the given
 * point under the map.
 * \param v  Second Cartesian coordinate  of the image of  the given
 * point under the map.
 *
 */
template <typename Mesh>
void PPS<Mesh>::from_canonical_domain(Halfedge *h, double x, double y,
                                      double &u, double &v) const {

#ifdef DEBUGMODE
  ERROR_UNLESS(h != 0, "Expected a non-null pointer");
#endif

  /*
   * Get the degree of the origin vertex of the given half-edge.
   */
  const unsigned int nu = get_degree(h);

  /*
   * Get the rotation angle associated with the given half-edge.
   */
  const double au = (2 * _MYPI) / nu;

  /*
   * Get the identifier of the given half-edge.
   */
  const unsigned i = get_id(h);

  /*
   * Change from Cartesian to polar coordinates.
   */
  double rr = std::sqrt((x * x) + (y * y));

  if (std::abs(rr) > 1e-16) {
    double aa;

    if (x < 0) {
      if (y >= 0) {
        aa = _MYPI - std::acos(std::abs(x) / rr);
      } else {
        aa = (2 * _MYPI) - std::acos(std::abs(x) / rr);
      }
    } else {
      if (y >= 0) {
        aa = std::asin(y / rr);
      } else {
        aa = -std::asin(std::abs(y) / rr);
      }
    }

    aa *= (6. / double(nu));

    rr *= (std::cos(_MYPI / nu) / std::cos(_MYPI / 6));

    u = rr * std::cos(aa);
    v = rr * std::sin(aa);
  } else {
    u = v = 0.;
  }

  /*
   * Perform  a 2D  rotation by  i  * au  around the  origin of  the
   * p-domain  associated  with  the  origin  vertex  of  the  given
   * half-edge.
   */
  rot_2d(u, v, (i * au), u, v);

  if (std::abs(u) < 1e-15)
    u = 0.;

  if (std::abs(v) < 1e-15)
    v = 0.;

  if (std::abs(1 - u) < 1e-15)
    u = 1.;

  if (std::abs(1 - v) < 1e-15)
    v = 1.;

  return;
}

/**
 * \fn void PPS< Mesh >::rot_2d( double x , double y , double ang , double& u ,
 * double& v ) const
 *
 * \brief Computes  the location of  a given point after  a rotation
 * around the origin by a given angle.
 *
 * \param x First Cartesian coordinate of the given point.
 * \param y Second Cartesian coordinate of the given point.
 * \param ang A given rotation angle.
 * \param u  First Cartesian coordinate of the  point resulting from
 * the rotation of (x,y) by ang around the origin.
 * \param v Second Cartesian  coordinate of the point resulting from
 * the rotation of (x,y) by ang around the origin.
 *
 */
template <typename Mesh>
void PPS<Mesh>::rot_2d(double x, double y, double ang, double &u,
                       double &v) const {
  u = (x * std::cos(ang)) - (y * std::sin(ang));
  v = (x * std::sin(ang)) + (y * std::cos(ang));
}

/**
 * \fn  unsigned int PPS< Mesh >::get_id( Halfedge* h ) const
 *
 * \brief  Returns  the  identifier   of  a  given  half-edge.   The
 * identifier of  a half-edge is a number  ranging from 0 to  N - 1,
 * where N is the degree of the origin vertex of the half-edge. Each
 * such a  number represents  the position of  one half-edge  in the
 * sequence   defined  by  a   counterclockwise  traversal   of  all
 * half-edges with  origin at the  same vertex. The  first half-edge
 * visited  in this  traversal  gets the  identifier  0, the  second
 * half-edge gets the identifier 1, and so on so forth.
 *
 * \param h Pointer to a half-edge.
 *
 * \returns The identifier of the given half-edge.
 *
 */
template <typename Mesh> unsigned int PPS<Mesh>::get_id(Halfedge *h) const {

#ifdef DEBUGMODE
  ERROR_UNLESS(h != 0, "Expected a non-null pointer");
#endif

  unsigned int i = 0;
  Halfedge *h2 = get_halfedge(get_org(h));
  while (h2 != h) {
    ++i;
    h2 = get_mate(get_prev(h2));
  }

  return i;
}

/**
 * \fn double PPS< Mesh >::weight_function( double x , double y , double R )
 * const
 *
 * \brief Evaluates the weight function at a given point.
 *
 * \param x First Cartesian coordinate of the given point.
 * \param y Second Cartesian coordinate of the given point.
 * \param R Radius of the circular support of the function.
 *
 * \return The function value at the given point.
 */
template <typename Mesh>
double PPS<Mesh>::weight_function(double x, double y, double R) const {
  double ll = std::sqrt((x * x) + (y * y));

  return eta_function(ll, 0.25 * R, R);
}

/**
 * \fn double PPS< Mesh >::eta_function( double s , double d1 , double d2 )
 * const
 *
 * \brief Evaluates  the eta blending function at  a given parameter
 * value.
 *
 * \param s Parameter value.
 * \param d1 Lower bound for the extent of the blending region.
 * \param d2 Upper bound for the extent of the blending region
 *
 * \return The function value at the given parameter value.
 */
template <typename Mesh>
double PPS<Mesh>::eta_function(double s, double d1, double d2) const {

#ifdef DEBUGMODE
  ERROR_UNLESS(d2 > d1, "Upper bound for the extent of blending region must be "
                        "larger than lower bound");
  ERROR_UNLESS(
      d1 > 0, "Lower bound for the extent of blending region must be positive");
  ERROR_UNLESS(
      d2 < 1,
      "Upper bound for the extent of blending region must be smaller than 1");
#endif

  double res;

  if (s <= d1) {
    res = 1;
  } else if (s < d2) {
    double h1 = (s - d1) / (d2 - d1);
    double h2 = 1 / std::sqrt(1 - h1);

    h1 = 1 / std::sqrt(h1);
    h2 = exp(h2 - h1);

    res = 1 / (1 + (h2 * h2));
  } else {
    res = 0;
  }

  return res;
}

} // namespace pps

/** @} */ // end of group class.
