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
#include <igl/write_triangle_mesh.h>
#include <modules/consistent_face_flippping.h>
#include <modules/remove_small_components.h>
#include <utils/mrf_potts.h>


// Mesh flipping.
DEFINE_bool(centerize, false, "");
DEFINE_bool(flip_x, false, "");
DEFINE_bool(flip_y, false, "");
DEFINE_bool(flip_z, false, "");

// Define input variables.
DEFINE_bool(run_face_labeling, false, "");
DEFINE_bool(run_part_disassembly, false, "");
DEFINE_bool(run_component_disassembly, false, "");
DEFINE_bool(run_point_sampling, false, "");
DEFINE_bool(run_barycenter_coloring, false, "");
DEFINE_bool(run_contacting_point_labeling, false, "");

// Face labeling params.
DEFINE_string(point_set, "", "point set file.");
DEFINE_string(point_labels, "", "point label file.");
DEFINE_string(out_mesh, "", "output mesh file.");
DEFINE_string(out_face_labels, "", "output face label file.");

// Part disassembly params.
DEFINE_string(out_part_mesh_dir, "", "output part mesh directory.");
DEFINE_string(out_part_mesh_unnormalized_dir, "",
    "output unnormalized part mesh directory.");
DEFINE_string(out_part_mesh_face_map_dir, "",
    "directory of output part mesh face map to input mesh.");

// Component disassembly params.
DEFINE_string(out_component_mesh_dir, "", "output component mesh directory.");
DEFINE_string(out_component_mesh_unnormalized_dir, "",
    "output unnormalized component mesh directory.");
//DEFINE_string(out_component_mesh_face_map_dir, "",
//    "directory of output component mesh face map to input mesh.");
DEFINE_string(out_component_mesh_face_map_file, "",
              "file of output component mesh face map to input mesh.");
DEFINE_int32(min_num_components, 3, "");
DEFINE_double(min_component_bbox_diagonal, 0.05, "");
DEFINE_bool(find_symmetric_components, false, "");

// Point sampling params.
DEFINE_int32(num_points, 1024, "");
DEFINE_string(out_point_set_dir, "", "output point set directory.");
DEFINE_bool(align_pca, false, "");
DEFINE_string(out_pca_alignment_dir, "",
    "directory of PCA alignment transformation matrix.");
DEFINE_bool(center_origin, false, "");
DEFINE_string(out_position_dir, "", "directory of center positions.");

// Barycenter-based mesh coloring params.
DEFINE_string(coloring_reference_mesh, "", "");

// Contacting point labeling params.
DEFINE_string(out_point_labels, "", "output point label file.");
DEFINE_double(max_contacting_squared_distance, 0.05*0.05, "");

// Transform mesh.
DEFINE_string(transformation, "", "");


bool LibiglMesh::read_point_set(const std::string& _filename, MatrixXd* _P) {
  CHECK_NOTNULL(_P);

  if (!Utils::read_eigen_matrix_from_file(_filename, _P, ' ')) {
    return false;
  }
  return true;
}

bool LibiglMesh::read_point_labels(
    const std::string& _filename, VectorXi* _PL) {
  CHECK_NOTNULL(_PL);

  if (!Utils::read_eigen_matrix_from_file(_filename, _PL)) {
    return false;
  }
  return true;
}

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

void LibiglMesh::processing() {
	bool mesh_modified = false;

	// Centerize.
	if (FLAGS_centerize) {
		V_ = V_.rowwise() - center_.transpose();
		mesh_modified = true;
	}

	// Mesh flipping.
	if (FLAGS_flip_x) { V_.col(0) = -V_.col(0); mesh_modified = true; }
	if (FLAGS_flip_y) { V_.col(1) = -V_.col(1); mesh_modified = true; }
	if (FLAGS_flip_z) { V_.col(2) = -V_.col(2); mesh_modified = true; }

	if (mesh_modified) {
		update_bounding_box();
    renderer_->set_mesh(V_, F_);
    renderer_->set_scene_pos(center_.cast<float>(), (float)radius_);
	}

  if (FLAGS_run_face_labeling) {
    processing_subdivide_mesh(FLAGS_out_mesh);

    processing_project_pts_labels_to_mesh(
        FLAGS_point_set,
        FLAGS_point_labels,
        FLAGS_out_face_labels);
  }
  else if (FLAGS_run_part_disassembly) {
    processing_disassemble_to_parts(
        FLAGS_out_part_mesh_dir,
        FLAGS_out_part_mesh_unnormalized_dir,
        FLAGS_out_part_mesh_face_map_dir);
  }
  else if (FLAGS_run_component_disassembly) {
    processing_disassemble_to_components(
        FLAGS_out_component_mesh_dir,
        FLAGS_out_component_mesh_unnormalized_dir,
        FLAGS_out_component_mesh_face_map_file,
        FLAGS_min_num_components,
        FLAGS_min_component_bbox_diagonal,
        FLAGS_find_symmetric_components);
  }
  else if (FLAGS_run_point_sampling) {
    processing_sample_points(
        FLAGS_num_points,
        FLAGS_out_point_set_dir,
        FLAGS_out_pca_alignment_dir,
        FLAGS_out_position_dir,
        FLAGS_align_pca,
        FLAGS_center_origin);
  }
  else if (FLAGS_run_barycenter_coloring) {
    processing_color_barycenter(FLAGS_coloring_reference_mesh);
  }
  else if (FLAGS_run_contacting_point_labeling) {
    processing_label_contacting_points(
        FLAGS_point_set,
        FLAGS_out_point_labels,
        FLAGS_max_contacting_squared_distance);
  }
  else if (FLAGS_transformation != "") {
    transform_mesh(FLAGS_transformation);
  }
}
