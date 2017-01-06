// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LibiglMeshViewer.h"

#include <igl/png/writePNG.h>


LibiglMeshViewer::LibiglMeshViewer(const int _width, const int _height)
  : LibiglMeshRendererT(_width, _height) {
  // FIXME(mhsung): Set window size.
}

LibiglMeshViewer::~LibiglMeshViewer() {
}

const Matrix4f& LibiglMeshViewer::get_projection() const {
  return viewer_.core.proj;
}

const Matrix4f& LibiglMeshViewer::get_modelview() const {
  return viewer_.core.model;
}

void LibiglMeshViewer::set_projection(const Matrix4f& _projection) {
  viewer_.core.proj = _projection;
}

void LibiglMeshViewer::set_modelview(const Matrix4f& _modelview) {
  viewer_.core.model = _modelview;
}

void LibiglMeshViewer::set_mesh(
    const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F) {
  viewer_.data.set_mesh(_V, _F);
}

void LibiglMeshViewer::set_face_colors(const Eigen::MatrixXf& _FC) {
  viewer_.data.set_colors(_FC.cast<double>());
}

void LibiglMeshViewer::run_loop() {
  viewer_.launch();
}

bool LibiglMeshViewer::snapshot(const std::string& _filename) {
  // Draw mesh.
  viewer_.draw();

  // Allocate temporary buffers.
  Matrix<unsigned char, Dynamic, Dynamic> R(width_, height_);
  Matrix<unsigned char, Dynamic, Dynamic> G(width_, height_);
  Matrix<unsigned char, Dynamic, Dynamic> B(width_, height_);
  Matrix<unsigned char, Dynamic, Dynamic> A(width_, height_);

  // Draw the scene in the buffers
  viewer_.core.draw_buffer(viewer_.data, viewer_.opengl, false, R, G, B, A);

  // Save it to a PNG.
  return igl::png::writePNG(R, G, B, A, _filename + std::string(".png"));
}
