// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LibiglMesh.h"

#include <Eigen/Geometry>
#include <igl/random_points_on_mesh.h>
#include <modules/PCA.h>
#include <utils/utils.h>


void LibiglMesh::sample_points_on_mesh(const int num_points) {
  SparseMatrix<double> B;
  VectorXi FI;
  igl::random_points_on_mesh(num_points, V_, F_, B, FI);
  P_ = B * V_;
}

void LibiglMesh::centerize_points(const std::string& _out_file) {
  const int num_samples = P_.rows();
  CHECK_GT(num_samples, 0);

  // Compute center and bounding box diagonal.
  const RowVector3d center = P_.colwise().mean();
  const auto bb_min = P_.colwise().minCoeff();
  const auto bb_max = P_.colwise().maxCoeff();
  const double bbox_diagonal = (bb_max - bb_min).norm();

  // Centerize point set.
  P_ = P_.rowwise() - center;

  if (_out_file != "") {
    RowVector4d center_and_size;
    center_and_size << center, bbox_diagonal;
    Utils::write_eigen_matrix_to_file(_out_file, center_and_size);
  }
}

void LibiglMesh::pca_align_points(const std::string& _out_file) {
  const int num_samples = P_.rows();
  CHECK_GT(num_samples, 0);

  // Compute PCA transformation matrix.
  Affine3d T = Affine3d::Identity();
  igl::PCA(P_, T);

  // PCA-align point set.
  const Eigen::Matrix<double, Dynamic, 3>& P_temp = P_;
  P_ = (T * P_temp.transpose()).transpose();

  const AngleAxisd R(T.rotation());
  double angle = R.angle();
  while (angle < 0) angle += 2 * M_PI;
  while (angle > 2 * M_PI) angle -= 2 * M_PI;
  Vector3d axis = R.axis().normalized();
  const Vector3d r = angle * axis;

  const Vector3d t = T.translation();

  // Scale along the first PCA axis.
  const double s = P_.col(0).maxCoeff() - P_.col(0).minCoeff();

  RowVectorXd transformation(7);
  transformation << r, t, s;

  if (_out_file != "") {
    Utils::write_eigen_matrix_to_file(_out_file, transformation);
  }
}
