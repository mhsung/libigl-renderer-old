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
#include <igl/point_mesh_squared_distance.h>
#include <igl/write_triangle_mesh.h>
#include <modules/consistent_face_flippping.h>
#include <modules/edge_lengths_simple.h>
#include <modules/PCA.h>
#include <modules/remove_duplicates_custom.h>
#include <modules/remove_small_components.h>
#include <modules/symmetric_elements.h>
#include <modules/upsample_non_manifold.h>
#include <utils/mrf_potts.h>


// Last modified: 01-09-2017
void LibiglMesh::processing_project_points_labels_to_mesh(
    const std::string& _out_face_labels_file) {
  const int num_points = P_.rows();
  const VectorXi label_set = Utils::unique(PL_);
  const int num_labels = label_set.size();

  // Construct a map from labels to label indices.
  std::unordered_map<int, int> label_idxs;
  for (int i = 0; i < num_labels; ++i) {
    label_idxs[ label_set[i] ] = i;
  }

  // Sample points on mesh.
  const int kAvgNumSamplesPerFace = 3;
  SparseMatrix<double> B;
  VectorXi FI;
  igl::random_points_on_mesh(kAvgNumSamplesPerFace * n_faces(), V_, F_, B, FI);
  const MatrixXd S = B * V_;
  const int num_samples = S.rows();
  LOG(INFO) << "# samples: " << num_samples;
  CHECK_GT(num_samples, 0);

  // Assign labels based on closest points.
  const KdTree::KdTreeXd P_tree(P_);
  VectorXi S_to_P_idxs;
  VectorXd S_to_P_sqr_dists;
  KdTree::find_k_closest_points(P_tree, S, &S_to_P_idxs, &S_to_P_sqr_dists);
  const VectorXi SL = Utils::slice_rows(PL_, S_to_P_idxs);

  // Compute CRF unary terms.
  MatrixXi F_to_S_counts = MatrixXi::Zero(n_faces(), num_labels);
  for (int i = 0; i < num_samples; ++i) {
    const int fid = FI[i];
    const int label_index = label_idxs[ SL[i] ];
    F_to_S_counts(fid, label_index) += 1;
  }

  // Normalize CRF unary terms and weight based on face areas.
  MatrixXd CRFUnary = MatrixXd::Zero(n_faces(), num_labels);
  VectorXd FA;
  igl::doublearea(V_, F_, FA);
  const VectorXd FA_normalized = FA.normalized();
  for (int fid = 0; fid < n_faces(); ++fid) {
    if (F_to_S_counts.row(fid).sum() > 0) {
      // Compute probability with Laplacian smoothing.
      RowVectorXd probs = (F_to_S_counts.row(fid).array() + 1).cast<double>();
      probs.normalize();
      CRFUnary.row(fid) = -probs.array().log();
    }
    CRFUnary.row(fid) *= FA_normalized[fid];
  }
  LOG(INFO) << "Max. Unary energy: " << CRFUnary.maxCoeff();

  // Compute binary penalties based on edge lengths.
  //igl::per_face_normals(V_, F_, FN_);
  //const double kCosFaceAngleTol = std::cos(60.0 / 180.0 * M_PI);
  const double kBinaryPenaltyWeight = 10.0;
  MatrixXi EV;
  MatrixXi FE;
  std::vector<std::list<int> > EF;
  VectorXd EL;
  igl::edge_topology_non_manifold(F_, EV, FE, EF);
  igl::edge_lengths_simple(V_, EV, EL);
  std::vector<Eigen::Triplet<double> > edge_coeffs;
  double max_penalty = 0.0;
  for (int eid = 0; eid < EF.size(); ++eid) {
    for (const int fid1 : EF[eid]) {
      for (const int fid2 : EF[eid]) {
        if (fid1 < fid2) {
          //if (FN_.row(fid1).dot(FN_.row(fid2)) < kCosFaceAngleTol) {
          //  continue;
          //}
          const double penalty = kBinaryPenaltyWeight * EL[eid];
          edge_coeffs.emplace_back(fid1, fid2, penalty);
          max_penalty = std::max(penalty, max_penalty);
        }
      }
    }
  }
  LOG(INFO) << "Max. Binary energy: " << max_penalty;

  /*
  // Enforce symmetric faces to have the same label.
  // Symmetry axis: Z-axis.
  const double kMRFMax = 1.0E6;
  Matrix3d ST = Matrix3d::Identity();
  ST(2, 2) = -1.0;
  VectorXi SFI;
  igl::symmetric_faces(V_, F_, ST, SFI);
  for (int fid = 0; fid < n_faces(); ++fid) {
    if (SFI[fid] > fid) {
      edge_coeffs.emplace_back(fid, SFI[fid], kMRFMax);
    }
  }
  */

  // Compute CRF binary terms.
  SparseMatrix<double> CRFBinary(n_faces(), n_faces());
  CRFBinary.setFromTriplets(edge_coeffs.begin(), edge_coeffs.end());

  // Run MRF solver and set face labels.
  const int kMaxMRFIters = 100;
  const VectorXi face_label_idxs = MRFEigen::solve_mrf_potts(
      CRFUnary, CRFBinary, MRFEigen::BP, kMaxMRFIters);
  FL_ = Utils::slice_rows(label_set, face_label_idxs);

  set_face_label_colors();

  if (_out_face_labels_file != "") {
    write_face_labels(_out_face_labels_file);
  }
}


void LibiglMesh::processing_MRF_with_point_label_probs(
        const std::string& _in_point_label_probs_file,
        const std::string& _out_face_labels_file) {
  const int num_points = P_.rows();
  MatrixXd P_unary_probs;
  if (!Utils::read_eigen_matrix_from_file(
        _in_point_label_probs_file, &P_unary_probs, ' ')) {
    LOG(ERROR) << "Can't read the file: '" <<
      _in_point_label_probs_file << "'";
    return;
  }
  CHECK_EQ(P_unary_probs.rows(), num_points);
  const int num_labels = P_unary_probs.cols();

  MatrixXd F_unary_probs = MatrixXd::Zero(n_faces(), num_labels);
  for (int fid = 0; fid < n_faces(); ++fid) {
    VectorXd sqrD;
    VectorXi I;
    MatrixXd C;
    igl::point_mesh_squared_distance(P_, V_, F_.row(fid), sqrD, I, C);
    const double kSquaredTol = std::pow(0.005, 2);

    int count = 0;
    for (int pid = 0; pid < num_points; ++pid) {
      if (sqrD[pid] < kSquaredTol) {
        F_unary_probs.row(fid) += P_unary_probs.row(pid);
        ++count;
      }
    }
    if (count > 0) F_unary_probs.row(fid) /= double(count);
  }

  // Set minimum probability.
  const double kMinProb = 1.0e-6;
  F_unary_probs = F_unary_probs.cwiseMax(kMinProb);

  // Compute unary penalties.
  MatrixXd CRFUnary = -F_unary_probs.array().log();

  // Weight based on face areas.
  VectorXd FA;
  igl::doublearea(V_, F_, FA);
  const VectorXd FA_normalized = FA.normalized();
  for (int fid = 0; fid < n_faces(); ++fid) {
    CRFUnary.row(fid) *= FA_normalized[fid];
  }

  // Compute binary penalties.
  igl::per_face_normals(V_, F_, FN_);
  const double kBinaryWeight = 1.0;
  const double kAngleWeight = 5.0;
  const double kMaxAngle = M_PI / 2.0;

  MatrixXi EV;
  MatrixXi FE;
  std::vector<std::list<int> > EF;
  igl::edge_topology_non_manifold(F_, EV, FE, EF);
  std::vector<Eigen::Triplet<double> > edge_coeffs;

  for (int eid = 0; eid < EF.size(); ++eid) {
    for (const int fid1 : EF[eid]) {
      for (const int fid2 : EF[eid]) {
        if (fid1 < fid2) {
          const auto cos_angle = FN_.row(fid1).dot(FN_.row(fid2));
          const auto angle = std::acos(
              std::min(std::max(cos_angle, -1.0), +1.0));
          const auto normalized_angle = std::min(angle, kMaxAngle) / kMaxAngle;
          const auto penalty = kBinaryWeight *
            std::exp(-kAngleWeight * normalized_angle);
          edge_coeffs.emplace_back(fid1, fid2, penalty);
        }
      }
    }
  }

  // Compute CRF binary terms.
  SparseMatrix<double> CRFBinary(n_faces(), n_faces());
  CRFBinary.setFromTriplets(edge_coeffs.begin(), edge_coeffs.end());

  // Run MRF solver and set face labels.
  const int kMaxMRFIters = 100;
  FL_ = MRFEigen::solve_mrf_potts(
      CRFUnary, CRFBinary, MRFEigen::TRW_S, kMaxMRFIters);

  /*
  FL_ = Eigen::VectorXi::Zero(n_faces());
  for(int fid = 0; fid < n_faces(); ++fid) {
    F_unary_probs.row(fid).maxCoeff(&FL_[fid]);
  }
  */

  // NOTE:
  // Do not render points.
  renderer_->set_points(Eigen::MatrixXd::Zero(0, 3));

  FL_ = FL_.array() + 1;
  set_face_label_colors();

  if (_out_face_labels_file != "") {
    write_face_labels(_out_face_labels_file);
  }
}
