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
    const int _num_points,
    const std::string& _out_point_set_dir,
    const std::string& _out_pca_alignment_dir,
    const std::string& _out_position_dir,
    const bool _align_pca,
    const bool _center_origin) {
  // Sample points on mesh.
  SparseMatrix<double> B;
  VectorXi FI;
  igl::random_points_on_mesh(_num_points, V_, F_, B, FI);
  Matrix<double, Dynamic, 3> P = B * V_;
  const int num_samples = P.rows();
  LOG(INFO) << "# samples: " << num_samples;
  CHECK_GT(num_samples, 0);

  Affine3d T = Affine3d::Identity();
	RowVector3d center = P.colwise().mean();
  double bbox_diagonal = 0.0;
  if (_align_pca) {
    igl::PCA(P, T);
    P = (T * P.transpose()).transpose();
  }
	else if (_center_origin) {
    P = P.rowwise() - center;
    const auto bb_min = P.colwise().minCoeff();
    const auto bb_max = P.colwise().maxCoeff();
    bbox_diagonal = (bb_max - bb_min).norm();
  }

  if (_out_point_set_dir != "") {
    const std::string kPointSetExt = ".pts";
    const std::string point_set_filename =
        filesystem::path(mesh_name_).basename() + kPointSetExt;
    const filesystem::path point_set_file =
        filesystem::path(_out_point_set_dir) /
            filesystem::path(point_set_filename);
    Utils::write_eigen_matrix_to_file(point_set_file.str(), P, ' ');
  }

  if (_align_pca && _out_pca_alignment_dir != "") {
    const std::string kMatrixExt = ".csv";
    const std::string transformation_filename =
        filesystem::path(mesh_name_).basename() + kMatrixExt;
    const filesystem::path transformation_file =
        filesystem::path(_out_pca_alignment_dir) /
        filesystem::path(transformation_filename);
    Utils::write_eigen_matrix_to_file(transformation_file.str(), T.matrix());
  }
	else if (_center_origin && _out_position_dir != "") {
    const std::string kMatrixExt = ".csv";
    const std::string center_filename =
        filesystem::path(mesh_name_).basename() + kMatrixExt;
    const filesystem::path center_file =
        filesystem::path(_out_position_dir) /
        filesystem::path(center_filename);
    RowVector4d center_and_size;
    center_and_size << center, bbox_diagonal;
    Utils::write_eigen_matrix_to_file(center_file.str(), center_and_size);
  }
}
