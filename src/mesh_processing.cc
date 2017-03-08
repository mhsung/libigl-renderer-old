// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LibiglMesh.h"

#include <Eigen/Geometry>
#include <utils/utils.h>


void LibiglMesh::normalize_mesh(MatrixXd& _V) {
  const auto bb_min = _V.colwise().minCoeff();
  const auto bb_max = _V.colwise().maxCoeff();
  const auto center = 0.5 * (bb_max + bb_min);
  const double bbox_diagonal = (bb_max - bb_min).norm();
  CHECK_GT(bbox_diagonal, 1.0E-6);

  // Move center to (0,0,0).
  _V = _V.rowwise() - center;

  // Scale to bounding box diagonal 1.
  _V /= bbox_diagonal;
}

void LibiglMesh::transform_mesh(const std::string& _filename) {
  Matrix4d mat;
  if (!Utils::read_eigen_matrix_from_file(_filename, &mat)) {
    return;
  }

  Affine3d T(mat);
  Matrix<double, Dynamic, 3> V_copy = V_;
  V_ = (T * V_copy.transpose()).transpose();

  update_bounding_box();
  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_mesh(V_, F_);
    renderer_->set_scene_pos(center_.cast<float>(), (float)radius_);
  }
}
