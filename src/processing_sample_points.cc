// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LibiglMesh.h"

#include <list>
#include <unordered_map>
#include <igl/barycenter.h>
#include <igl/doublearea.h>
#include <igl/facet_components.h>
#include <igl/random_points_on_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <igl/per_face_normals.h>
#include <igl/write_triangle_mesh.h>
#include <modules/consistent_face_flippping.h>
#include <modules/edge_lengths_simple.h>
#include <modules/PCA.h>
#include <modules/remove_duplicates_custom.h>
#include <modules/remove_small_components.h>
#include <modules/symmetric_elements.h>
#include <modules/upsample_non_manifold.h>
#include <utils/filesystem/path.h>
#include <utils/mrf_potts.h>


void LibiglMesh::processing_sample_points(
    const int _num_sample_points,
    const std::string& _out_point_set_dir,
    const std::string& _out_pca_transformation_dir,
    const std::string& _out_position_dir,
    const std::string& _out_centered_point_set_dir,
    const std::string& _out_pca_aligned_point_set_dir) {
  // Sample points on mesh.
  SparseMatrix<double> B;
  VectorXi FI;
  igl::random_points_on_mesh(_num_sample_points, V_, F_, B, FI);
  const Matrix<double, Dynamic, 3> P = B * V_;
  const int num_samples = P.rows();
  LOG(INFO) << "# samples: " << num_samples;
  CHECK_GT(num_samples, 0);

  // Compute center and bounding box diagonal.
	const RowVector3d center = P.colwise().mean();
  const auto bb_min = P.colwise().minCoeff();
  const auto bb_max = P.colwise().maxCoeff();
  const double bbox_diagonal = (bb_max - bb_min).norm();

  // Compute PCA transformation matrix.
  Affine3d T = Affine3d::Identity();
  igl::PCA(P, T);


  const std::string basename = filesystem::path(mesh_name_).basename();

  if (_out_point_set_dir != "") {
    const filesystem::path out_file =
      filesystem::path(_out_point_set_dir) /
      filesystem::path(basename + std::string(".pts"));

    Utils::write_eigen_matrix_to_file(out_file.str(), P, ' ');
  }

	if (_out_position_dir != "") {
    const filesystem::path out_file =
      filesystem::path(_out_position_dir) /
      filesystem::path(basename + std::string(".csv"));

    RowVector4d center_and_size;
    center_and_size << center, bbox_diagonal;
    Utils::write_eigen_matrix_to_file(out_file.str(), center_and_size);
  }

  if (_out_pca_transformation_dir != "") {
    const filesystem::path out_file =
      filesystem::path(_out_pca_transformation_dir) /
      filesystem::path(basename + std::string(".csv"));

    Utils::write_eigen_matrix_to_file(out_file.str(), T.matrix());
  }

  if (_out_centered_point_set_dir != "") {
    const filesystem::path out_file =
      filesystem::path(_out_centered_point_set_dir) /
      filesystem::path(basename + std::string(".pts"));

    const Matrix<double, Dynamic, 3> centered_P = P.rowwise() - center;
    Utils::write_eigen_matrix_to_file(out_file.str(), centered_P, ' ');
  }

  if (_out_pca_aligned_point_set_dir != "") {
    const filesystem::path out_file =
      filesystem::path(_out_pca_aligned_point_set_dir) /
      filesystem::path(basename + std::string(".pts"));

    const Matrix<double, Dynamic, 3> aligned_P =
      (T * P.transpose()).transpose();
    Utils::write_eigen_matrix_to_file(out_file.str(), aligned_P, ' ');
  }
}
