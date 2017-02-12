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


void LibiglMesh::processing_disassemble_to_components(
    const std::string& _out_component_mesh_dir,
    const std::string& _out_component_mesh_unnormalized_dir,
    const std::string& _out_component_mesh_face_map_dir,
    const int _min_num_components,
    const int _min_component_bbox_diagonal) {
  // Find components.
  igl::facet_components(F_, FL_);

  const VectorXi label_set = Utils::unique(FL_);
  if (label_set.size() < _min_num_components) {
    LOG(WARNING) << "Warning: Too few components. Skip.";
    return;
  }

  // Create normalized output mesh directory.
  const filesystem::path component_mesh_dir(_out_component_mesh_dir);
  if (_out_component_mesh_dir != "" &&
      !component_mesh_dir.is_directory()) {
    CHECK(filesystem::create_directory(component_mesh_dir));
  }

  // Create unnormalized output mesh directory.
  const filesystem::path component_mesh_unnormalized_dir(
      _out_component_mesh_unnormalized_dir);
  if (_out_component_mesh_unnormalized_dir != "" &&
      !component_mesh_unnormalized_dir.is_directory()) {
    CHECK(filesystem::create_directory(component_mesh_unnormalized_dir));
  }

  // Create output face map directory.
  const filesystem::path component_mesh_face_map_dir(
      _out_component_mesh_face_map_dir);
  if (_out_component_mesh_face_map_dir != "" &&
      !component_mesh_face_map_dir.is_directory()) {
    CHECK(filesystem::create_directory(component_mesh_face_map_dir));
  }

  // Process each component.
  const int kNumDigits = 4;
  const std::string kMeshExt = ".obj";
  const std::string kFaceMapExt = ".txt";
  int label = 0;

  for (int i = 0; i < label_set.size(); ++i) {
    const VectorXi label_fids = Utils::find(FL_, label_set[i]);
    MatrixXi label_F_old = Utils::slice_rows(F_, label_fids);
    MatrixXd label_V;
    MatrixXi label_F;
    VectorXi IX;
    igl::remove_unreferenced(V_, label_F_old, label_V, label_F, IX);

    // Ignore if the component is too small.
    const auto bb_min = label_V.colwise().minCoeff();
    const auto bb_max = label_V.colwise().maxCoeff();
    const double bbox_diagonal = (bb_max - bb_min).norm();
    LOG(INFO) << "Label (" << label << "): BBox diagonal = " << bbox_diagonal;
    if (bbox_diagonal < _min_component_bbox_diagonal) continue;

    std::stringstream label_mesh_sstr;
    label_mesh_sstr << std::setw(kNumDigits) << std::setfill('0') << label;
    ++label;

    // Write unnormalized component mesh.
    if (_out_component_mesh_unnormalized_dir != "") {
      const filesystem::path component_mesh_file =
          filesystem::path(_out_component_mesh_unnormalized_dir) /
          filesystem::path(label_mesh_sstr.str() + kMeshExt);
      igl::write_triangle_mesh(component_mesh_file.str(), label_V, label_F);
    }

    // Write normalize component mesh.
    normalize_mesh(label_V);

    if (_out_component_mesh_dir != "") {
      const filesystem::path component_mesh_file =
          filesystem::path(_out_component_mesh_dir) /
          filesystem::path(label_mesh_sstr.str() + kMeshExt);
      igl::write_triangle_mesh(component_mesh_file.str(), label_V, label_F);
    }

    // Write component mesh face map to input mesh.
    if (_out_component_mesh_face_map_dir != "") {
      const filesystem::path component_mesh_face_map_file =
          filesystem::path(_out_component_mesh_face_map_dir) /
          filesystem::path(label_mesh_sstr.str() + kFaceMapExt);
      Utils::write_eigen_matrix_to_file(
          component_mesh_face_map_file.str(), label_fids);
    }
  }
}
