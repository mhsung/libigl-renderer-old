// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LibiglMeshT.h"

#include <igl/read_triangle_mesh.h>
#include <utils/filesystem/path.h>
#include <utils/google_tools.h>
#include <utils/utils.h>

// Define input variables.
DEFINE_string(mesh, "", "mesh file.");
DEFINE_string(face_labels, "", "face label file.");
DEFINE_double(azimuth_deg, 0.0, "azimuth (degree). "
    "ignored if 'modelview_matrix' is set");
DEFINE_double(elevation_deg, 0.0, "elevation (degree). "
    "ignored if 'modelview_matrix' is set");
DEFINE_double(theta_deg, 0.0, "theta (degree). "
    "ignored if 'modelview_matrix' is set");
DEFINE_string(projection_matrix, "", "projection matrix file.");
DEFINE_string(modelview_matrix, "", "modelview matrix file.");
DEFINE_string(bbox, "", "bounding box file.");
DEFINE_string(snapshot, "", "snapshot file.");


LibiglMeshT::LibiglMeshT()
  : renderer_(nullptr) {
}

LibiglMeshT::LibiglMeshT(LibiglMeshRendererT* _renderer)
  : renderer_(_renderer) {
}

bool LibiglMeshT::read_mesh(const std::string& _filename) {
  if (!igl::read_triangle_mesh(_filename, V_, F_)) {
    LOG(WARNING) << "Can't read the file: '" << _filename << "'";
    return false;
  }
  mesh_name_ = filesystem::path(_filename).filename();

  update_bounding_box();

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_mesh(V_, F_);
    renderer_->set_scene_pos(center_.cast<float>(), (float)radius_);
  }

  return true;
}

bool LibiglMeshT::read_face_labels(const std::string& _filename) {
  FL_ = VectorXi(n_faces());
  if (!Utils::read_eigen_matrix_from_file(_filename, &FL_)) {
    return false;
  }

  // Set face colors.
  set_face_label_colors();
  return true;
}

bool LibiglMeshT::write_face_labels(const std::string& _filename) {
  if (!Utils::write_eigen_matrix_to_file(_filename, FL_)) {
    return false;
  }
  return true;
}

void LibiglMeshT::set_face_label_colors() {
  if (FL_.rows() != F_.rows()) {
    LOG(WARNING) << "Number of face labels does not match number of faces.";
    return;
  }

  FC_ = MatrixXf(n_faces(), 3);
  for (int fid = 0; fid < n_faces(); ++fid) {
    Vector3f color;
    Utils::random_label_rgb_color(FL_(fid), &color);
    FC_.row(fid) = color.transpose();
  }

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_face_colors(FC_);
  }
}

void LibiglMeshT::update_bounding_box() {
  bb_min_ = V_.colwise().minCoeff();
  bb_max_ = V_.colwise().maxCoeff();
  center_ = 0.5 * (bb_min_ + bb_max_);
  radius_ = 0.5 * (bb_max_ - bb_min_).norm();
}

bool LibiglMeshT::write_bounding_box(const std::string& _filename) {
  const double bbox_diagonal = (bb_max_ - bb_min_).norm();
	Eigen::VectorXd bb_info(4);
  bb_info << center_, bbox_diagonal;
  if (!Utils::write_eigen_matrix_to_file(_filename, bb_info.transpose())) {
    return false;
  }
  return true;
}

void LibiglMeshT::pre_processing() {
  if (!read_mesh(FLAGS_mesh)) {
    return;
  }

  if (FLAGS_face_labels != "") {
    if (!read_face_labels(FLAGS_face_labels)) {
      return;
    }
  }
}

void LibiglMeshT::post_processing() {
  if (FLAGS_azimuth_deg != 0.0 || FLAGS_elevation_deg != 0.0 ||
      FLAGS_theta_deg != 0.0) {
    if (renderer_ == nullptr) {
      LOG(WARNING) << "Renderer is not set";
    } else {
      renderer_->set_camera_params(
          Vector3f(FLAGS_azimuth_deg, FLAGS_elevation_deg, FLAGS_theta_deg),
          (float)radius_);
    }
  }

  if (FLAGS_projection_matrix != "") {
    if (renderer_ == nullptr) {
      LOG(WARNING) << "Renderer is not set";
    } else if (!renderer_->read_projection(FLAGS_projection_matrix)) {
      return;
    }
  }

  if (FLAGS_modelview_matrix != "") {
    if (renderer_ == nullptr) {
      LOG(WARNING) << "Renderer is not set";
    } else if (!renderer_->read_modelview(FLAGS_modelview_matrix)) {
      return;
    }
  }

	if (FLAGS_bbox != "") {
		if (!write_bounding_box(FLAGS_bbox)) {
			return;
		}
	}

#ifdef USE_OSMESA
  renderer_->snapshot(FLAGS_snapshot);
#else
  renderer_->run_loop();
#endif
}

void LibiglMeshT::parse_arguments_and_run() {
  pre_processing();
  processing();
  post_processing();
}
