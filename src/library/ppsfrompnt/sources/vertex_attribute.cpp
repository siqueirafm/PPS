/**
 * \file vertex_attribute.cpp
 *
 * \brief Implementation of class \c VertexAttribute.
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

#include "vertex_attribute.hpp"

namespace ppsfrompnt {

using pps::Bezier;

VertexAttribute::VertexAttribute() : _patch(nullptr) {}

VertexAttribute::VertexAttribute(Bezier *patch) : _patch(patch) {}

VertexAttribute::VertexAttribute(const VertexAttribute &a) {
  _patch = new Bezier(*(a.get_patch()));
}

VertexAttribute::VertexAttribute(VertexAttribute &&a) : _patch(a._patch) {
  a._patch = nullptr;
}

VertexAttribute::~VertexAttribute() {
  if (get_patch())
    delete get_patch();
}

VertexAttribute &VertexAttribute::operator=(const VertexAttribute &a) {
  if (get_patch())
    delete get_patch();

  _patch = new Bezier(*(a.get_patch()));

  return *this;
}

Bezier *VertexAttribute::get_patch() const { return _patch; }

void VertexAttribute::set_patch(Bezier *patch) { _patch = patch; }

} // namespace ppsfrompnt
