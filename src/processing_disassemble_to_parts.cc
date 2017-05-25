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


void LibiglMesh::processing_disassemble_to_parts(
    const std::string& _out_part_mesh_dir,
    const std::string& _out_part_mesh_unnormalized_dir,
    const std::string& _out_part_mesh_face_map_dir) {
  // Check whether face labels are loaded.
  CHECK_EQ(F_.rows(), FL_.rows());

  const VectorXi label_set = Utils::unique(FL_);
  const int num_digits = std::to_string(label_set.maxCoeff()).size();

  // Create root output mesh directory.
  const filesystem::path part_mesh_dir(_out_part_mesh_dir);
  if (_out_part_mesh_dir != "" && !part_mesh_dir.is_directory()) {
    CHECK(filesystem::create_directory(part_mesh_dir));
  }

  // Create root output mesh directory.
  const filesystem::path part_mesh_unnormalized_dir(
      _out_part_mesh_unnormalized_dir);
  if (_out_part_mesh_unnormalized_dir != "" &&
      !part_mesh_unnormalized_dir.is_directory()) {
    CHECK(filesystem::create_directory(part_mesh_unnormalized_dir));
  }

  // Create root output face map directory.
  const filesystem::path part_mesh_face_map_dir(
      _out_part_mesh_face_map_dir);
  if (_out_part_mesh_face_map_dir != "" &&
      !part_mesh_face_map_dir.is_directory()) {
    CHECK(filesystem::create_directory(part_mesh_face_map_dir));
  }

  // Process each part.
  for (int i = 0; i < label_set.size(); ++i) {
    const int label = label_set[i];
    const VectorXi label_fids = Utils::find(FL_, label);
    MatrixXi label_F_old = Utils::slice_rows(F_, label_fids);
    MatrixXd label_V;
    MatrixXi label_F;
    VectorXi IX;
    igl::remove_unreferenced(V_, label_F_old, label_V, label_F, IX);

    // Remove small components.
    igl::remove_small_components(label_V, label_F);

    std::stringstream label_sstr;
    label_sstr << std::setw(num_digits) << std::setfill('0') << label;

    // Write unnormalized part mesh.
    if (_out_part_mesh_unnormalized_dir != "") {
      const filesystem::path part_mesh_file = part_mesh_unnormalized_dir /
        filesystem::path(label_sstr.str() + std::string(".obj"));
      igl::write_triangle_mesh(part_mesh_file.str(), label_V, label_F);
    }

    // Write normalized part mesh.
    normalize_mesh(label_V);

    if (_out_part_mesh_dir != "") {
      const filesystem::path part_mesh_file = part_mesh_dir /
        filesystem::path(label_sstr.str() + std::string(".obj"));
      igl::write_triangle_mesh(part_mesh_file.str(), label_V, label_F);
    }

    // Write part mesh face map to input mesh.
    if (_out_part_mesh_face_map_dir != "") {
      const filesystem::path part_mesh_face_map_file = part_mesh_face_map_dir /
        filesystem::path(label_sstr.str() + std::string(".seg"));
      Utils::write_eigen_matrix_to_file(
          part_mesh_face_map_file.str(), label_fids);
    }
  }
}
