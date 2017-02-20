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
#include <igl/hausdorff.h>
#include <igl/random_points_on_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <igl/per_face_normals.h>
#include <igl/write_triangle_mesh.h>
#include <modules/consistent_face_flippping.h>
#include <modules/edge_lengths_simple.h>
#include <modules/merge_meshes.h>
#include <modules/PCA.h>
#include <modules/remove_duplicates_custom.h>
#include <modules/remove_small_components.h>
#include <modules/symmetric_elements.h>
#include <modules/upsample_non_manifold.h>
#include <utils/filesystem/path.h>
#include <utils/mrf_potts.h>
#include <utils/SparseICP/ICP.h>


const int kNumSegmentPointSamples = 1000;
const double kMinIOU = 0.50;
const double kMinAlignmentError = 5.0;


struct MeshSegment {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  MatrixXd V_;
  MatrixXi F_;
  MatrixXd FN_;
  VectorXi fids_;
  MatrixXd P_aligned_;
  MatrixXd PN_aligned_;

  void merge(const MeshSegment& _other) {
    igl::add_mesh(V_, F_, _other.V_, _other.F_);
    const VectorXi fids_copy = fids_;
    fids_.resize(fids_copy.size() + _other.fids_.size());
    fids_ << fids_copy, _other.fids_;
  }

  void update_face_normals() {
    igl::per_face_normals(V_, F_, Vector3d(1,1,1).normalized(), FN_);
  }

  void update_aligned_point_samples(const int _num_points) {
    SparseMatrix<double> B;
    VectorXi FI;
    igl::random_points_on_mesh(_num_points, V_, F_, B, FI);
    const Matrix<double, Dynamic, 3> P = B * V_;

    Affine3d T = Affine3d::Identity();
    igl::PCA(P, T);
    P_aligned_ = (T * P.transpose()).transpose();

    // Compute point normals.
    update_face_normals();
    const MatrixXd PN = Utils::slice_rows(FN_, FI);
    PN_aligned_ = (T.rotation() * PN.transpose()).transpose();
  }
};

typedef std::unique_ptr<MeshSegment> MeshSegmentPtr;


double point_alignment_error(const Matrix3Xd& P1, const Matrix3Xd& P2) {
  nanoflann::KDTreeAdaptor<MatrixBase<Matrix<double, 3, -1, 0, 3, -1>>, 3,
      nanoflann::metric_L2_Simple> kdtree_2(P2);

  double max_dist = 0;
  for(int i = 0; i < P1.cols(); ++i) {
    // Find the closest points from P1 to P2.
    int idx12;
    double sq_dist12;
    kdtree_2.query(P1.col(i).data(), 1, &idx12, &sq_dist12);

    // Compute neighbor distance at the point of P2.
    const int kNumNeighbors = std::min(10, (int)P2.cols());
    VectorXi idx22(kNumNeighbors);
    VectorXd sq_dist22(kNumNeighbors);
    kdtree_2.query(P2.col(idx12).data(), kNumNeighbors, idx22.data(),
                   sq_dist22.data());

    // Compute alignment error over neighbor distance.
    const double dist = std::sqrt(sq_dist12) / std::sqrt(sq_dist22.maxCoeff());
    max_dist = std::max(max_dist, dist);
  }

  return max_dist;
}

bool is_segment_symmetric(
    const MeshSegment& _mesh_1, const MeshSegment& _mesh_2) {
  if (_mesh_1.V_.rows() != _mesh_2.V_.rows() ||
      _mesh_1.F_.rows() != _mesh_2.F_.rows()) {
    return false;
  }
  /*
  // Compare PCA-aligned bounding box.
  const RowVectorXd bb_min_1 = _mesh_1.P_aligned_.colwise().minCoeff();
  const RowVectorXd bb_max_1 = _mesh_1.P_aligned_.colwise().maxCoeff();
  const RowVectorXd bb_size_1 = bb_max_1 - bb_min_1;

  const RowVectorXd bb_min_2 = _mesh_2.P_aligned_.colwise().minCoeff();
  const RowVectorXd bb_max_2 = _mesh_2.P_aligned_.colwise().maxCoeff();
  const RowVectorXd bb_size_2 = bb_max_2 - bb_min_2;

  const RowVectorXd min_bb_size = bb_size_1.cwiseMin(bb_size_2);
  const RowVectorXd max_bb_size = bb_size_1.cwiseMax(bb_size_2);
  const double intersection_volume = min_bb_size.prod();
  const double union_volume = max_bb_size.prod();
  const double iou = (union_volume > 1.0E-6) ?
                     intersection_volume / union_volume : 0.0;
  if (iou < kMinIOU) return false;

  // Compute Hausdorff distance.
  // P1 and P2 will be modified.
  Matrix3Xd P1 = _mesh_1.P_aligned_.transpose();
  Matrix3Xd P2 = _mesh_2.P_aligned_.transpose();
  ICP::point_to_point(P1, P2);
  const double alignment_error = std::max(
      point_alignment_error(P1, P2), point_alignment_error(P2, P1));
  if (alignment_error > kMinAlignmentError) return false;
   */

  return true;
}

std::vector<MeshSegmentPtr> merge_symmetric_segments(
    const std::vector<MeshSegmentPtr>& _segments) {
  const int num_segments = _segments.size();

  // Initial cluster IDs of segments.
  VectorXi segment_cids;
  igl::colon(0, 1, num_segments - 1, segment_cids);

  for (int i = 0; i < num_segments - 1; ++i) {
    for (int j = i + 1; j < num_segments; ++j) {
      if (is_segment_symmetric(*_segments[i], *_segments[j])) {
        // Merge clusters.
        const int min_cid_ij = std::min(segment_cids(i), segment_cids(j));
        for (int k = 0; k < num_segments; ++k) {
          if (segment_cids(k) == segment_cids(i) ||
              segment_cids(k) == segment_cids(j)) {
            segment_cids(k) = min_cid_ij;
          }
        }
      }
    }
  }

  const VectorXi cid_set = Utils::unique(segment_cids);
  std::vector<MeshSegmentPtr> merged_segments;
  merged_segments.reserve(cid_set.size());

  for (int i = 0; i < cid_set.size(); ++i) {
    MeshSegmentPtr cluster_mesh(new MeshSegment());

    // Get segment IDs of i-th cluster.
    const VectorXi cluster_sids = Utils::find(segment_cids, cid_set[i]);
    for (int j = 0; j < cluster_sids.size(); ++j) {
      const int sid = cluster_sids[j];
      cluster_mesh->merge(*_segments[sid]);
    }

    merged_segments.push_back(std::move(cluster_mesh));
  }

  return merged_segments;
}

void LibiglMesh::processing_disassemble_to_components(
    const std::string& _out_component_mesh_dir,
    const std::string& _out_component_mesh_unnormalized_dir,
    const std::string& _out_component_mesh_face_map_file,
    const int _min_num_components,
    const int _min_component_bbox_diagonal,
    const bool _find_symmetry) {
  // Find components.
  igl::facet_components(F_, FL_);

  const VectorXi label_set = Utils::unique(FL_);
  if (label_set.size() < _min_num_components) {
    LOG(WARNING) << "Warning: Too few components. Skip.";
    return;
  }

  // Process each component.
  std::vector<MeshSegmentPtr> label_meshes;
  label_meshes.reserve(label_set.size());

  for (int i = 0; i < label_set.size(); ++i) {
    MeshSegmentPtr label_mesh(new MeshSegment());

    label_mesh->fids_ = Utils::find(FL_, label_set[i]);
    MatrixXi label_F_old = Utils::slice_rows(F_, label_mesh->fids_);
    VectorXi IX;
    igl::remove_unreferenced(
        V_, label_F_old, label_mesh->V_, label_mesh->F_, IX);

    // Ignore if the component is too small.
    const auto bb_min = label_mesh->V_.colwise().minCoeff();
    const auto bb_max = label_mesh->V_.colwise().maxCoeff();
    const double bbox_diagonal = (bb_max - bb_min).norm();
    if (bbox_diagonal < _min_component_bbox_diagonal) continue;

    label_meshes.push_back(std::move(label_mesh));
  }

  // Update segment point samples.
//  for (auto& label_mesh : label_meshes) {
//    label_mesh->update_aligned_point_samples(kNumSegmentPointSamples);
//  }

  // Merge symmetric components.
  if (_find_symmetry) {
    std::vector<MeshSegmentPtr> symmetric_label_meshes =
        merge_symmetric_segments(label_meshes);
    label_meshes.swap(symmetric_label_meshes);
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

  /*
  // Create output face map directory.
  const filesystem::path component_mesh_face_map_dir(
      _out_component_mesh_face_map_dir);
  if (_out_component_mesh_face_map_dir != "" &&
      !component_mesh_face_map_dir.is_directory()) {
    CHECK(filesystem::create_directory(component_mesh_face_map_dir));
  }
  */

  const int kNumDigits = 4;
  const std::string kMeshExt = ".obj";
  const std::string kFaceMapExt = ".txt";

  // NOTE:
  // 'Zero' label indicates unassigned.
  int label = 1;
  for (auto& label_mesh : label_meshes) {
    std::stringstream label_mesh_sstr;
    label_mesh_sstr << std::setw(kNumDigits) << std::setfill('0') << label;
    ++label;

    // Assign component IDs.
    for (int i = 0; i < label_mesh->fids_.size(); ++i) {
      const int fid = label_mesh->fids_[i];
      FL_(fid) = label;
    }

    // Write unnormalized component mesh.
    if (_out_component_mesh_unnormalized_dir != "") {
      const filesystem::path component_mesh_file =
          filesystem::path(_out_component_mesh_unnormalized_dir) /
          filesystem::path(label_mesh_sstr.str() + kMeshExt);
      igl::write_triangle_mesh(component_mesh_file.str(),
                               label_mesh->V_, label_mesh->F_);
    }

    // Write normalize component mesh.
    normalize_mesh(label_mesh->V_);

    if (_out_component_mesh_dir != "") {
      const filesystem::path component_mesh_file =
          filesystem::path(_out_component_mesh_dir) /
          filesystem::path(label_mesh_sstr.str() + kMeshExt);
      igl::write_triangle_mesh(component_mesh_file.str(),
                               label_mesh->V_, label_mesh->F_);
    }

    /*
    // Write component mesh face map to input mesh.
    if (_out_component_mesh_face_map_dir != "") {
      const filesystem::path component_mesh_face_map_file =
          filesystem::path(_out_component_mesh_face_map_dir) /
          filesystem::path(label_mesh_sstr.str() + kFaceMapExt);
      Utils::write_eigen_matrix_to_file(
          component_mesh_face_map_file.str(), label_mesh->fids_);
    }
    */
  }

  // Write component ID face labels.
  if (_out_component_mesh_face_map_file != "") {
    write_face_labels(_out_component_mesh_face_map_file);
  }

  set_face_label_colors();
}
