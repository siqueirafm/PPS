/**
 * \file face_attribute.cpp
 *
 * \brief Implementation of class \c FaceAttribute.
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

#include "face_attribute.hpp"

namespace ppsfrompnt {

FaceAttribute::FaceAttribute() : _patch(nullptr) {}

FaceAttribute::FaceAttribute(PNTriangle *patch) : _patch(patch) {}

FaceAttribute::FaceAttribute(const FaceAttribute &a) {
  _patch = new PNTriangle(*(a.get_patch()));
}

FaceAttribute::FaceAttribute(FaceAttribute &&a) : _patch(a._patch) {
  a._patch = nullptr;
}

FaceAttribute::~FaceAttribute() {
  if (get_patch())
    delete get_patch();
}

FaceAttribute &FaceAttribute::operator=(const FaceAttribute &a) {
  if (get_patch())
    delete get_patch();

  _patch = new PNTriangle(*(a.get_patch()));

  return *this;
}

PNTriangle *FaceAttribute::get_patch() const { return _patch; }

void FaceAttribute::set_patch(PNTriangle *patch) { _patch = patch; }

} // namespace ppsfrompnt
