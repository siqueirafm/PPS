#pragma once

/**
 * \file halfedge_attribute.hpp
 *
 * \brief Definition of  the class HalfedgeAttribute, which represents
 * a  set of halfedge  attributes of  the underlying  triangle surface
 * mesh  of a Parametric  Pseudo-Surface (PPS)  constructed from  a PN
 * triangle surface.
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

#include <halfedge.hpp> // dcel::Halfedge

#include "face_attribute.hpp"   // FaceAttribute
#include "vertex_attribute.hpp" // VertexAttribute

/**
 * \defgroup PPSFROMPNTNameSpace Namespace ppsfrompnt.
 * @{
 */

/**
 * \namespace ppsfrompnt
 *
 * \brief  The  namespace  ppsfrompnt   contains  the  definition  and
 * implementation   of   all   classes   representing   a   Parametric
 * Pseudo-Surface (PPS) that  approximates the geometry of  a given PN
 * triangle surface.
 */

namespace ppsfrompnt {

/**
 * \class HalfedgeAttribute
 *
 * \brief This class  represents a set of attributes  of a half-edge
 * of the underlying triangle surface mesh of a PPS constructed from
 * a PN triangle surface.
 */
class HalfedgeAttribute {
public:
  /**
   * \using Halfedge
   *
   * \brief  Defines  Halfedge   as  an  alias  for  dcel::Halfedge<
   * VertexAttribute, FaceAttribute , int , HalfedgeAttribute >
   */
  using Halfedge =
      dcel::Halfedge<VertexAttribute, FaceAttribute, int, HalfedgeAttribute>;

  /**
   * \fn HalfedgeAttribute()
   *
   * \brief Creates an instance of this class.
   */
  HalfedgeAttribute();

  /**
   * \fn HalfedgeAttribute( Halfedge* h )
   *
   * \brief Creates an instance of this class.
   *
   * \param h A pointer to the halfedge that owns this attribute.
   */
  explicit HalfedgeAttribute(Halfedge *h);

  /**
   * \fn Halfedge* get_owner() const
   *
   * \brief  Returns  a pointer  to  the  half-edge  that owns  this
   * attribute.
   *
   * \return A pointer to the half-edge that owns this attribute.
   */
  Halfedge *get_owner() const;

  /**
   * \fn void set_owner( Halfedge* h )
   *
   * \brief Assigns an address to  the pointer to the half-edge that
   * owns this attribute.
   *
   * \param h A pointer to a half-edge.
   */
  void set_owner(Halfedge *h);

  /**
   * \fn unsigned get_pps_id()
   *
   * \brief Returns the PPS identifier (ID) of this half-edge.  This
   * identifier is a number from 0  to N-1, where N is the degree of
   * the origin vertex of this halfedge. The number assigned to this
   * half-edge is  its position  in a counterclockwise  traversal of
   * all half-edges incident to the origin vertex of this half-edge.
   * The  first  half-edge  of   this  traversal  is  the  half-edge
   * associated with the vertex and its ID is equal to 0. The second
   * half-edge gets the ID 1, and  so on. This method is a modifier,
   * as the  ID is computed  and stored in  a data member  the first
   * time the method is invoked.
   *
   * \return The PPS identifier of this half-edge.
   */
  unsigned get_pps_id();

  /**
   * \fn void set_pps_id( unsigned id )
   *
   * \brief Assigns a value to  the data member representing the PPS
   * identifier of this half-edge.
   *
   * \param  id An  unsigned  integer representing  a half-edge  PPS
   * identifier.
   */
  void set_pps_id(unsigned id);

  /**
   * \fn bool get_pps_id_flag() const
   *
   * \brief Returns  the logic value  true if the PPS  identifier of
   * this  half-edge has already  been computed;  otherwise, returns
   * the logic value false.
   *
   * \return  The logic  value true  if the  PPS identifier  of this
   * half-edge  has already  been computed;  otherwise,  returns the
   * logic value false.
   */
  bool get_pps_id_flag() const;

  /**
   * \fn void set_pps_id_flag( bool flag )
   *
   * \brief Assigns a  logic value (true or false)  to the flag that
   * indicates if  the PPS identifier of this  half-edge has already
   * been computed.
   *
   * \param flag A  logic value ( true or false )  to be assigned to
   * the PPS identifier flag of this half-edge.
   */
  void set_pps_id_flag(bool flag);

  /**
   * \fn unsigned get_origin_vertex_degree()
   *
   * \brief  Returns  the  degree  of  the  origin  vertex  of  this
   * half-edge. If  the degree has  not been computed  yet, computes
   * and  returns the  degree, which  means  that this  method is  a
   * modifier.
   *
   * \return The degree of this vertex.
   */
  unsigned get_origin_vertex_degree();

  /**
   * \fn void set_origin_vertex_degree( unsigned degree )
   *
   * \brief  Assigns a  value to  the data  member  representing the
   * degree of the origin vertex of this half-edge.
   *
   * \param degree An unsigned integer representing a vertex degree.
   */
  void set_origin_vertex_degree(unsigned degree);

  /**
   * \fn bool get_degree_flag() const
   *
   * \brief Returns the logic value true if the degree of the origin
   * vertex of this half-edge  has already been computed; otherwise,
   * returns the logic value false.
   *
   * \return The logic value true if the degree of the origin vertex
   * of this half-edge has already been computed; otherwise, returns
   * the logic value false.
   */
  bool get_degree_flag() const;

  /**
   * \fn void set_degree_flag( bool flag )
   *
   * \brief Assigns a  logic value (true or false)  to the flag that
   * indicates if the degree of  the origin vertex of this half-edge
   * has already been computed.
   *
   * \param flag A  logic value ( true or false )  to be assigned to
   * the degree flag of this half-edge.
   */
  void set_degree_flag(bool flag);

private:
  /**
   * \fn unsigned compute_pps_id() const
   *
   * \brief  Computes and returns  the PPS  identifier (ID)  of this
   * half-edge.
   *
   * \return The PPS identifier of this half-edge.
   */
  unsigned compute_pps_id() const;

  /**
   * \fn unsigned compute_origin_vertex_degree() const
   *
   * \brief Computes and returns the  degree of the origin vertex of
   * this half-edge.
   *
   * \return The degree of the origin vertex of this half-edge.
   */
  unsigned compute_origin_vertex_degree() const;

private:
  Halfedge *_owner; ///< Pointer to the half-edge that owns this attribute..
  unsigned _id;  ///< PPS identifier of the half-edge that owns this attribute.
  bool _id_flag; ///< PPS identifier flag of the half-edge that owns this
                 ///< attribute..
  unsigned _degree;  ///< Degree of the origin vertex of the half-edge that owns
                     ///< this attribute.
  bool _degree_flag; ///< Degree flag of the half-edge that owns this attribute.
};

} // namespace ppsfrompnt

/** @} */ // end of group class.
