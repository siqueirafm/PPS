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
 */

#include "face_attribute.hpp"

namespace ppsfromloop {

FaceAttribute::FaceAttribute() : _patch(nullptr) {}

FaceAttribute::FaceAttribute(LOOPatch *patch) : _patch(patch) {}

FaceAttribute::FaceAttribute(const FaceAttribute &a) {
  _patch = new LOOPatch(*(a.get_patch()));
}

FaceAttribute::FaceAttribute(FaceAttribute &&a) : _patch(a._patch) {
  a._patch = nullptr;
}

FaceAttribute::~FaceAttribute() {
  if (_patch)
    delete _patch;
}

FaceAttribute &FaceAttribute::operator=(const FaceAttribute &a) {
  if (get_patch())
    delete get_patch();

  _patch = new LOOPatch(*(a.get_patch()));

  return *this;
}

LOOPatch *FaceAttribute::get_patch() const { return _patch; }

void FaceAttribute::set_patch(LOOPatch *patch) { _patch = patch; }

} // namespace ppsfromloop
