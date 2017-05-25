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
#include <utils/mrf_potts.h>


// Last modified: 01-09-2017
void LibiglMesh::processing_project_pts_labels_to_mesh(
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
