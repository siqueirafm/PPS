#pragma once

/**
 * \file ppssampler.hpp
 *
 * \brief Definition  and implementation of a class  that represents a
 * Parametric  Pseudo-Surface  (PPS) sampler.  This  sampler yields  a
 * triangle  mesh  that approximates  the  PPS,  which  is useful  for
 * visualization purposes.
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

#include "pps.hpp" // pps::PPS

#include <exception.hpp> // ERROR_UNLESS

#include <map>    // std::map
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
 * \class PPSsampler
 *
 * \brief  This   class  represents  a  sampler   for  a  Parametric
 * Pseudo-Surface (PPS). This sampler  provides us with a method for
 * creating a  triangle mesh that approximates the  PPS. This method
 * can be used for visualization purposes.
 *
 */
template <typename Mesh> class PPSsampler {
public:
  // ---------------------------------------------------------------
  //
  // Type definitions
  //
  // ---------------------------------------------------------------

  /**
   * \using Vertex
   *
   * \brief Defines an alias for Mesh::Vertex.
   */
  using Vertex = typename Mesh::Vertex;

  /**
   * \using Halfedge
   *
   * \brief Defines an alias for Mesh::Halfedge.
   */
  using Halfedge = typename Mesh::Halfedge;

  /**
   * \using Face
   *
   * \brief Defines an alias for Mesh::Face.
   */
  using Face = typename Mesh::Face;

  /**
   * \using VertexIterator
   *
   * \brief Defines an alias for Mesh::VertexIterator.
   */
  using VertexIterator = typename Mesh::VertexIterator;

  /**
   * \using FaceIterator
   *
   * \brief Defines an alias for Mesh::FaceIterator.
   */
  using FaceIterator = typename Mesh::FaceIterator;

  // ---------------------------------------------------------------
  //
  // Public methods
  //
  // ---------------------------------------------------------------

  /**
   * \fn PPSsampler( PPS< Mesh >* pps )
   *
   * \brief Creates an instance of this class.
   *
   * \param pps A pointer to a PPS.
   *
   */
  PPSsampler(PPS<Mesh> *pps) : _pps(pps) {}

  /**
   * \fn virtual ~PPSsampler()
   *
   * \brief Virtual destructor.
   */
  virtual ~PPSsampler() {}

  /**
   * \fn void sample( unsigned lod , unsigned& numpts , double*& lvpps ,
   * double*& lvmsh , unsigned& numfaces , unsigned*& lfaces )
   *
   * \brief Sample this PPS  with a given Level-Of-Detail (LOD). The
   * sample  is  carried  out  in a  triangle  midpoint  subdivision
   * manner.
   *
   * \param lod The level-of-detail of the requested sample.
   * \param numpts The number of sample points.
   * \param lvpps A set of sample points on the PPS.
   * \param lvmsh A set of sample points on the surface approximated
   * by the PPS.
   * \param numfaces The number of faces.
   * \param lfaces A set of faces connecting the sample points.
   */
  void sample(unsigned lod, unsigned &numpts, double *&lvpps, double *&lvmsh,
              unsigned &numfaces, unsigned *&lfaces);

private:
  // ---------------------------------------------------------------
  //
  // Private methods
  //
  // ---------------------------------------------------------------

  /**
   * \fn void midpointsub( unsigned lod , unsigned& vcount , Face* face )
   *
   * \brief Samples a PPS and the surface approximated by it.
   *
   * \param lod Level-of-detail of the sampling.
   * \param vcount Vertex counter.
   * \param face Pointer to a face of the PPS.
   */
  void midpointsub(unsigned lod, unsigned &vcount, Face *face);

  /**
   * \fn void midpointsub( unsigned lod , unsigned& vcount , Face* face ,
   * unsigned v1 , unsigned v2 , unsigned v3 , double c1[ 3 ] , double c2[ 3 ] ,
   * double c3[ 3 ] )
   *
   * \brief Samples a PPS and the surface approximated by it.
   *
   * \param lod Level-of-detail of the sampling.
   * \param vcount Vertex counter.
   * \param face Pointer to a face of the PPS.
   * \param v1 ID of the first vertex of the face.
   * \param v2 ID of the second vertex of the face.
   * \param v3 ID of the third vertex of the face.
   * \param c1 Barycentric coordinate of the first vertex of the face.
   * \param c2 Barycentric coordinate of the second vertex of the face.
   * \param c3: barycentric coordinate of the third vertex of the face.
   */
  void midpointsub(unsigned lod, unsigned &vcount, Face *face, unsigned v1,
                   unsigned v2, unsigned v3, double c1[3], double c2[3],
                   double c3[3]);

  // ---------------------------------------------------------------
  //
  // Private data members
  //
  // ---------------------------------------------------------------

  pps::PPS<Mesh> *_pps; ///< A pointer to a PPS.

  std::vector<double> _lv1; ///< An array of sample points on the PPS.

  std::vector<double> _lv2; ///< An array of sample points on the surface
                            ///< approximated by the PPS.

  std::vector<unsigned>
      _lf; ///< An array of faces connecting the sample points.

  std::map<Vertex *, unsigned> _vnt; ///< A hash table for storing vertex ID's

  std::map<std::pair<unsigned, unsigned>, unsigned>
      _vpa; ///< A hash table for storing edge ID's
};

// -----------------------------------------------------------------
//
// Implementation of public methods
//
// -----------------------------------------------------------------

/**
 * \fn void PPSsampler< Mesh >::sample( unsigned lod , unsigned& numpts ,
 * double*& lvpps , double*& lvmsh , unsigned& numfaces , unsigned*& lfaces )
 *
 * \brief Sample  this PPS with  a given Level-Of-Detail  (LOD). The
 * sample is carried out in a triangle midpoint subdivision manner.
 *
 * \param lod The level-of-detail of the requested sample.
 * \param numpts The number of sample points.
 * \param lvpps A set of sample points on the PPS.
 * \param lvmsh A  set of sample points on  the surface approximated
 * by the PPS.
 * \param numfaces The number of faces.
 * \param lfaces A set of faces connecting the sample points.
 */
template <typename Mesh>
void PPSsampler<Mesh>::sample(unsigned lod, unsigned &numpts, double *&lvpps,
                              double *&lvmsh, unsigned &numfaces,
                              unsigned *&lfaces) {
  unsigned i = 0;

  for (VertexIterator vit = _pps->vertices_begin(); !_pps->is_done(vit);
       _pps->move_forward(vit)) {

    Vertex *vertex = *vit;

    _vnt.insert(std::make_pair(vertex, i));

    /**
     * Sample the PPS and the surface approximated by it.
     */

    double pt1[3];
    double pt2[3];

    Halfedge *he = vertex->get_halfedge();

    Face *face = he->get_face();

    if (he == face->get_halfedge()) {
      _pps->eval_pps(face, 1, 0, 0, pt1[0], pt1[1], pt1[2]);
      _pps->eval_surface(face, 1, 0, 0, pt2[0], pt2[1], pt2[2]);
    } else if (he->get_next() == face->get_halfedge()) {
      _pps->eval_pps(face, 0, 0, 1, pt1[0], pt1[1], pt1[2]);
      _pps->eval_surface(face, 0, 0, 1, pt2[0], pt2[1], pt2[2]);
    } else {
      _pps->eval_pps(face, 0, 1, 0, pt1[0], pt1[1], pt1[2]);
      _pps->eval_surface(face, 0, 1, 0, pt2[0], pt2[1], pt2[2]);
    }

    /**
     * Store the vertex in the set of vertices.
     */
    _lv1.push_back(pt1[0]);
    _lv1.push_back(pt1[1]);
    _lv1.push_back(pt1[2]);

    // Store the vertex in the set of vertices.
    _lv2.push_back(pt2[0]);
    _lv2.push_back(pt2[1]);
    _lv2.push_back(pt2[2]);

    ++i;
  }

  /**
   * Sample each  face using  a "midpoint" subdivision  (like Loop's
   * scheme), and then connect the points to form the faces. We keep
   * the points and faces in a list.
   */
  for (FaceIterator fit = _pps->faces_begin(); !_pps->is_done(fit);
       _pps->move_forward(fit)) {
    midpointsub(lod, i, *fit);
  }

  ERROR_UNLESS(
      _lv1.size() == _lv2.size(),
      "Number of vertices of loop and pps simplicial surfaces must agree");

  numpts = _lv1.size();
  numfaces = _lf.size();

  lvpps = new double[numpts];
  lvmsh = new double[numpts];
  lfaces = new unsigned[numfaces];

  for (i = 0; i < numpts; i++) {
    lvpps[i] = _lv1[i];
    lvmsh[i] = _lv2[i];
  }

  for (i = 0; i < numfaces; i++)
    lfaces[i] = _lf[i];

  _lv1.clear();
  _lv2.clear();

  _lf.clear();

  _vnt.clear();
  _vpa.clear();

  return;
}

/**
 * \fn void PPSsampler< Mesh >::midpointsub( unsigned lod , unsigned& vcount ,
 * Face* face )
 *
 * \brief Samples a PPS and the surface approximated by it.
 *
 * \param lod Level-of-detail of the sampling.
 * \param vcount Vertex counter.
 * \param face Pointer to a face of the PPS.
 */
template <typename Mesh>
void PPSsampler<Mesh>::midpointsub(unsigned lod, unsigned &vcount, Face *face) {
  /**
   * Get the vertices of the given face.
   */
  Vertex *v1 = face->get_halfedge()->get_origin();
  Vertex *v2 = face->get_halfedge()->get_next()->get_origin();
  Vertex *v3 = face->get_halfedge()->get_prev()->get_origin();

  /**
   * Get the ID of the above vertices.
   */
  unsigned id1 = _vnt.find(v1)->second;
  unsigned id2 = _vnt.find(v2)->second;
  unsigned id3 = _vnt.find(v3)->second;

  /**
   * If  we already reached  the largest  LOD then  we can  create a
   * face.
   */
  if (lod == 0) {
    _lf.push_back(id1);
    _lf.push_back(id2);
    _lf.push_back(id3);
  } else {
    double c1[3] = {1, 0, 0};
    double c2[3] = {0, 1, 0};
    double c3[3] = {0, 0, 1};

    midpointsub(lod - 1, vcount, face, id1, id2, id3, c1, c2, c3);
  }

  return;
}

/**
 * \fn void PPSsampler< Mesh >::midpointsub( unsigned lod , unsigned& vcount ,
 * Face* face , unsigned v1 , unsigned v2 , unsigned v3 , double c1[ 3 ] ,
 * double c2[ 3 ] , double c3[ 3 ] )
 *
 * \brief Samples a PPS and the surface approximated by it.
 *
 * \param lod Level-of-detail of the sampling.
 * \param vcount Vertex counter.
 * \param face Pointer to a face of the PPS.
 * \param v1 ID of the first vertex of the face.
 * \param v2 ID of the second vertex of the face.
 * \param v3 ID of the third vertex of the face.
 * \param c1 Barycentric coordinate of the first vertex of the face.
 * \param c2 Barycentric coordinate of the second vertex of the face.
 * \param c3: barycentric coordinate of the third vertex of the face.
 */
template <typename Mesh>
void PPSsampler<Mesh>::midpointsub(unsigned lod, unsigned &vcount, Face *face,
                                   unsigned v1, unsigned v2, unsigned v3,
                                   double c1[3], double c2[3], double c3[3]) {
  /**
   * Compute the midpoint vertices if they do not exist already.
   */
  double c12[3];
  double c23[3];
  double c31[3];

  for (unsigned i = 0; i < 3; i++) {
    c12[i] = 0.5 * (c1[i] + c2[i]);
    c23[i] = 0.5 * (c2[i] + c3[i]);
    c31[i] = 0.5 * (c3[i] + c1[i]);
  }

  std::map<std::pair<unsigned, unsigned>, unsigned>::iterator pit =
      _vpa.find(std::make_pair(v1, v2));

  unsigned v12;
  if (pit != _vpa.end()) {
    v12 = pit->second;
  } else {
    double pt1[3];
    double pt2[3];

    _pps->eval_pps(face, c12[0], c12[1], c12[2], pt1[0], pt1[1], pt1[2]);
    _pps->eval_surface(face, c12[0], c12[1], c12[2], pt2[0], pt2[1], pt2[2]);

    _lv1.push_back(pt1[0]);
    _lv1.push_back(pt1[1]);
    _lv1.push_back(pt1[2]);

    _lv2.push_back(pt2[0]);
    _lv2.push_back(pt2[1]);
    _lv2.push_back(pt2[2]);

    v12 = vcount;
    _vpa.insert(std::make_pair(std::make_pair(v2, v1), v12));
    ++vcount;
  }

  pit = _vpa.find(std::make_pair(v2, v3));

  unsigned v23;
  if (pit != _vpa.end()) {
    v23 = pit->second;
  } else {
    double pt1[3];
    double pt2[3];

    _pps->eval_pps(face, c23[0], c23[1], c23[2], pt1[0], pt1[1], pt1[2]);
    _pps->eval_surface(face, c23[0], c23[1], c23[2], pt2[0], pt2[1], pt2[2]);

    _lv1.push_back(pt1[0]);
    _lv1.push_back(pt1[1]);
    _lv1.push_back(pt1[2]);

    _lv2.push_back(pt2[0]);
    _lv2.push_back(pt2[1]);
    _lv2.push_back(pt2[2]);

    v23 = vcount;
    _vpa.insert(std::make_pair(std::make_pair(v3, v2), v23));
    ++vcount;
  }

  pit = _vpa.find(std::make_pair(v3, v1));

  unsigned v31;
  if (pit != _vpa.end()) {
    v31 = pit->second;
  } else {
    double pt1[3];
    double pt2[3];

    _pps->eval_pps(face, c31[0], c31[1], c31[2], pt1[0], pt1[1], pt1[2]);
    _pps->eval_surface(face, c31[0], c31[1], c31[2], pt2[0], pt2[1], pt2[2]);

    _lv1.push_back(pt1[0]);
    _lv1.push_back(pt1[1]);
    _lv1.push_back(pt1[2]);

    _lv2.push_back(pt2[0]);
    _lv2.push_back(pt2[1]);
    _lv2.push_back(pt2[2]);

    v31 = vcount;
    _vpa.insert(std::make_pair(std::make_pair(v1, v3), v31));
    ++vcount;
  }

  if (lod == 0) {
    _lf.push_back(v1);
    _lf.push_back(v12);
    _lf.push_back(v31);

    _lf.push_back(v12);
    _lf.push_back(v2);
    _lf.push_back(v23);

    _lf.push_back(v31);
    _lf.push_back(v23);
    _lf.push_back(v3);

    _lf.push_back(v12);
    _lf.push_back(v23);
    _lf.push_back(v31);
  } else {
    midpointsub(lod - 1, vcount, face, v1, v12, v31, c1, c12, c31);
    midpointsub(lod - 1, vcount, face, v12, v2, v23, c12, c2, c23);
    midpointsub(lod - 1, vcount, face, v31, v23, v3, c31, c23, c3);
    midpointsub(lod - 1, vcount, face, v12, v23, v31, c12, c23, c31);
  }

  return;
}

} // namespace pps

/** @} */ // end of group class.
