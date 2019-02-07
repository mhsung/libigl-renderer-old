// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "LibiglMeshT.h"

#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <utils/filesystem/path.h>
#include <utils/google_tools.h>
#include <utils/utils.h>

// Define input variables.
DEFINE_string(mesh, "", "mesh file.");
DEFINE_string(vertex_values, "", "vertex value file.");
DEFINE_string(vertex_labels, "", "vertex label file.");
DEFINE_string(face_labels, "", "face label file.");
DEFINE_string(point_set, "", "point set file.");
DEFINE_string(point_labels, "", "point label file.");
DEFINE_string(point_values, "", "point value file.");
DEFINE_string(keypoints, "", "keypoint coordinates.");
DEFINE_string(keypoint_labels, "", "keypoint labels.");
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
DEFINE_string(out_mesh, "", "output mesh file.");
DEFINE_string(out_point_set, "", "output point set file.");
DEFINE_string(out_point_labels, "", "output point label file.");


LibiglMeshT::LibiglMeshT()
  : renderer_(nullptr),
    mesh_name_(""),
    V_(MatrixXd(0, 3)),
    F_(MatrixXi(0, 3)),
    VC_(MatrixXf(0, 3)),
    FC_(MatrixXf(0, 3)),
    VL_(VectorXi(0, 3)),
    FL_(VectorXi(0, 3)),
    VN_(MatrixXd(0, 3)),
    FN_(MatrixXd(0, 3)),
    P_(MatrixXd(0, 3)),
    PC_(MatrixXf(0, 3)),
    PL_(VectorXi(0, 3)),
    PN_(MatrixXd(0, 3)),
    bb_min_(Vector3d::Zero()),
    bb_max_(Vector3d::Zero()),
    center_(Vector3d::Zero()),
    radius_(1.0) {
}

LibiglMeshT::LibiglMeshT(LibiglMeshRendererT* _renderer)
  : renderer_(_renderer) {
}

bool LibiglMeshT::read_mesh(const std::string& _filename) {
  if (!igl::read_triangle_mesh(_filename, V_, F_)) {
    LOG(WARNING) << "Can't read the file: '" << _filename << "'";
    return false;
  }

  if (mesh_name_ == "" ) mesh_name_ = filesystem::path(_filename).filename();
  update_bounding_box();

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_mesh(V_, F_);
    renderer_->set_scene_pos(center_.cast<float>(), (float)radius_);
  }

  return true;
}

bool LibiglMeshT::read_vertex_labels(const std::string& _filename) {
  if (!Utils::read_eigen_matrix_from_file(_filename, &VL_)) {
    return false;
  }

  // Set vertex colors.
  set_vertex_label_colors();
  return true;
}

bool LibiglMeshT::read_vertex_values(const std::string& _filename) {
  VectorXf VV;
  if (!Utils::read_eigen_matrix_from_file(_filename, &VV)) {
    return false;
  }

  if (VV.size() != V_.rows()) {
    LOG(WARNING) << "Number of vertex values does not match number of vertices.";
    return false;
  }

  VC_ = compute_color_map(VV);

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_vertex_colors(VC_);
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

bool LibiglMeshT::read_point_set(const std::string& _filename) {
  if (!Utils::read_eigen_matrix_from_file(_filename, &P_, ' ')) {
    LOG(WARNING) << "Can't read the file: '" << _filename << "'";
    return false;
  }

  if (mesh_name_ == "" ) mesh_name_ = filesystem::path(_filename).filename();
  update_bounding_box();

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_points(P_);
    renderer_->set_scene_pos(center_.cast<float>(), (float)radius_);
  }

  return true;
}

bool LibiglMeshT::read_point_labels(const std::string& _filename) {
  if (!Utils::read_eigen_matrix_from_file(_filename, &PL_)) {
    return false;
  }

  // Set point colors.
  set_point_label_colors();
  return true;
}

bool LibiglMeshT::read_point_values(const std::string& _filename) {
  VectorXf PV;
  if (!Utils::read_eigen_matrix_from_file(_filename, &PV)) {
    return false;
  }

  if (PV.size() != P_.rows()) {
    LOG(WARNING) << "Number of point values does not match number of points.";
    return false;
  }

  PC_ = compute_color_map(PV);

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_point_colors(PC_);
  }

  return true;
}

bool LibiglMeshT::read_keypoints(const std::string& _str) {
  /*
  MatrixXd KP_;
  if (!Utils::read_eigen_matrix_from_file(_filename, &KP_, ' ')) {
    LOG(WARNING) << "Can't read the file: '" << _filename << "'";
    return false;
  }
  */
  std::vector<std::string> strs = Utils::split_string(_str);
  const int n_keypoints = int(strs.size()) / 3;
  CHECK_EQ(3 * n_keypoints, strs.size());

  KP_ = MatrixXd::Zero(n_keypoints, 3);
  for (int i = 0; i < n_keypoints; ++i) {
    for (int j = 0; j < 3; ++j) {
      KP_(i, j) = std::stod(strs[3*i+j]);
    }
  }

  KPC_ = MatrixXf(n_keypoints, 3);
  for (int i = 0; i < n_keypoints; ++i) {
    Vector3f color;
    Utils::random_label_rgb_color(i + 1, &color);
    KPC_.row(i) = color.transpose();
  }

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_keypoints(KP_);
    renderer_->set_keypoint_colors(KPC_);
  }

  return true;
}

bool LibiglMeshT::read_keypoint_labels(const std::string& _str) {
  std::vector<std::string> strs = Utils::split_string(_str);
  const int n_keypoints = KP_.rows();
  CHECK_EQ(strs.size(), n_keypoints);

  KPL_ = VectorXi::Zero(n_keypoints);
  for (int i = 0; i < n_keypoints; ++i) {
    KPL_(i) = std::stoi(strs[i]);
  }

  KPC_ = MatrixXf(n_keypoints, 3);
  for (int i = 0; i < n_keypoints; ++i) {
    Vector3f color;
    Utils::random_label_rgb_color(KPL_[i], &color);
    KPC_.row(i) = color.transpose();
  }

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_keypoints(KP_);
    renderer_->set_keypoint_colors(KPC_);
  }

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

MatrixXf LibiglMeshT::compute_color_map(const VectorXf& _values) {
  const int n_values = _values.size();
  MatrixXf colors = MatrixXf(n_values, 3);
  colors.setZero();

  const float vmin = _values.minCoeff();
  const float vmax = _values.maxCoeff();
  const float dv = vmax - vmin;

  if (dv > 1.0e-8) {
    for (int pid = 0; pid < n_values; ++pid) {
      Vector3f color = Vector3f::Ones();
      const float v = _values[pid];

      // https://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale
      if (v < (vmin + 0.25f * dv)) {
        color[0] = 0.0f;
        color[1] = 4.0f * (v - vmin) / dv;
      } else if (v < (vmin + 0.5f * dv)) {
        color[0] = 0.0f;
        color[2] = 1.0f + 4.0f * (vmin + 0.25f * dv - v) / dv;
      } else if (v < (vmin + 0.75f * dv)) {
        color[0] = 4.0f * (v - vmin - 0.5f * dv) / dv;
        color[2] = 0.0f;
      } else {
        color[1] = 1.0f + 4.0f * (vmin + 0.75f * dv - v) / dv;
        color[2] = 0.0f;
      }

      colors.row(pid) = color.transpose();
    }
  }

  return colors;
}

void LibiglMeshT::set_vertex_label_colors() {
  if (VL_.rows() != V_.rows()) {
    LOG(WARNING) << "Number of vertex labels does not match number of vertexs.";
    return;
  }

  VC_ = MatrixXf(n_vertices(), 3);
  for (int pid = 0; pid < n_vertices(); ++pid) {
    Vector3f color;
    Utils::random_label_rgb_color(VL_(pid), &color);
    VC_.row(pid) = color.transpose();
  }

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_vertex_colors(VC_);
  }
}

void LibiglMeshT::set_point_label_colors() {
  if (PL_.rows() != P_.rows()) {
    LOG(WARNING) << "Number of point labels does not match number of points.";
    return;
  }

  PC_ = MatrixXf(n_points(), 3);
  for (int pid = 0; pid < n_points(); ++pid) {
    Vector3f color;
    Utils::random_label_rgb_color(PL_(pid), &color);
    PC_.row(pid) = color.transpose();
  }

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_point_colors(PC_);
  }
}

void LibiglMeshT::update_bounding_box() {
  MatrixXd points = MatrixXd::Zero(V_.rows() + P_.rows(), 3);
  if (points.rows() == 0) return;

  points << V_, P_;
  bb_min_ = points.colwise().minCoeff();
  bb_max_ = points.colwise().maxCoeff();
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
  if (FLAGS_mesh != "") {
    if (!read_mesh(FLAGS_mesh)) {
      return;
    }
  }

  if (FLAGS_vertex_labels != "") {
    if (!read_vertex_labels(FLAGS_vertex_labels)) {
      return;
    }
  }

  if (FLAGS_vertex_values != "") {
    if (!read_vertex_values(FLAGS_vertex_values)) {
      return;
    }
  }

  if (FLAGS_face_labels != "") {
    if (!read_face_labels(FLAGS_face_labels)) {
      return;
    }
  }

  if (FLAGS_point_set != "") {
    if (!read_point_set(FLAGS_point_set)) {
      return;
    }
  }

  if (FLAGS_point_labels != "") {
    if (!read_point_labels(FLAGS_point_labels)) {
      return;
    }
  }

  if (FLAGS_point_values != "") {
    if (!read_point_values(FLAGS_point_values)) {
      return;
    }
  }

  if (FLAGS_keypoints != "") {
    if (!read_keypoints(FLAGS_keypoints)) {
      return;
    }
  }

  if (FLAGS_keypoint_labels != "") {
    if (!read_keypoint_labels(FLAGS_keypoint_labels)) {
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

  if (FLAGS_out_mesh != "") {
    igl::write_triangle_mesh(FLAGS_out_mesh, V_, F_);
  }

  if (FLAGS_out_point_set != "") {
    if (P_.rows() == PN_.rows()) {
      const int n_points = P_.rows();
      CHECK_EQ(P_.cols(), 3);
      CHECK_EQ(PN_.cols(), 3);
      MatrixXd P_and_PN(n_points, 6);
      P_and_PN.leftCols(3) = P_;
      P_and_PN.rightCols(3) = PN_;
      Utils::write_eigen_matrix_to_file(FLAGS_out_point_set, P_and_PN, ' ');
    } else {
      Utils::write_eigen_matrix_to_file(FLAGS_out_point_set, P_, ' ');
    }
  }

  if (FLAGS_out_point_labels != "") {
    CHECK_EQ(PL_.rows(), P_.rows());
    Utils::write_eigen_matrix_to_file(FLAGS_out_point_labels, PL_, ' ');
  }

  // Compute normals.
  //igl::per_face_normals(V_, F_, FN_);
  //igl::per_vertex_normals(V_, F_, FN_, VN_);

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
