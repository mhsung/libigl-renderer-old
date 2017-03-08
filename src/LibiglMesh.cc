// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LibiglMesh.h"
#include <utils/utils.h>


// Mesh processing.
DEFINE_bool(centerize, false, "");
DEFINE_bool(flip_x, false, "");
DEFINE_bool(flip_y, false, "");
DEFINE_bool(flip_z, false, "");
DEFINE_string(translation, "", "");
DEFINE_string(transformation, "", "");

// Point set processing.
DEFINE_bool(sample_points, false, "");
DEFINE_int32(num_sample_points, 1024, "");
DEFINE_bool(centerize_point_set, false, "");
DEFINE_string(out_point_set_center, "", "");
DEFINE_bool(pca_align_point_set, false, "");
DEFINE_string(out_pca_transformation, "", "");

// Additional processing.
DEFINE_bool(run_face_labeling, false, "");
DEFINE_bool(run_part_disassembly, false, "");
DEFINE_bool(run_component_disassembly, false, "");
DEFINE_bool(run_point_sampling, false, "");
DEFINE_bool(run_point_transformation, false, "");
DEFINE_bool(run_barycenter_coloring, false, "");
DEFINE_bool(run_contacting_point_labeling, false, "");

// Face labeling params.
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

// Barycenter-based mesh coloring params.
DEFINE_string(coloring_reference_mesh, "", "");

// Contacting point labeling params.
DEFINE_string(out_point_labels, "", "output point label file.");
DEFINE_double(max_contacting_squared_distance, 0.005 * 0.005, "");


void LibiglMesh::mesh_processing() {
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

  if (FLAGS_translation != "") {
    std::vector<std::string> strs = Utils::split_string(FLAGS_translation);
    RowVector3d t;
    t << std::stod(strs[0]), std::stod(strs[1]), std::stod(strs[2]);
    V_ = V_.rowwise() + t;
    mesh_modified = true;
  }

  if (FLAGS_transformation != "") {
    transform_mesh(FLAGS_transformation);
    mesh_modified = true;
  }

  if (mesh_modified) {
    update_bounding_box();
    renderer_->set_mesh(V_, F_);
    renderer_->set_scene_pos(center_.cast<float>(), (float)radius_);
  }
}

void LibiglMesh::point_set_processing() {
  if (FLAGS_sample_points) {
    sample_points_on_mesh(FLAGS_num_sample_points);
  }

  if (FLAGS_centerize_point_set) {
    centerize_points(FLAGS_out_point_set_center);
  }

  if (FLAGS_pca_align_point_set) {
    pca_align_points(FLAGS_out_pca_transformation);
  }
}

void LibiglMesh::processing() {
  mesh_processing();
  point_set_processing();

  if (FLAGS_run_face_labeling) {
    processing_subdivide_mesh();
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
  else if (FLAGS_run_barycenter_coloring) {
    processing_color_barycenter(FLAGS_coloring_reference_mesh);
  }
  else if (FLAGS_run_contacting_point_labeling) {
    processing_label_contacting_points(
        FLAGS_point_set,
        FLAGS_out_point_labels,
        FLAGS_max_contacting_squared_distance);
  }
}
