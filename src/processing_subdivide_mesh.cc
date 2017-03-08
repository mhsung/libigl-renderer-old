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
#include <modules/consistent_face_flippping.h>
#include <modules/edge_lengths_simple.h>
#include <modules/PCA.h>
#include <modules/remove_duplicates_custom.h>
#include <modules/remove_small_components.h>
#include <modules/symmetric_elements.h>
#include <modules/upsample_non_manifold.h>
#include <utils/mrf_potts.h>


void LibiglMesh::remove_duplicates() {
  LOG(INFO) << "Before removing duplicates:";
  LOG(INFO) << " - # vertices: " << n_vertices();
  LOG(INFO) << " - # faces: " << n_faces();

  igl::remove_duplicate_vertices_custom(V_, F_);
  LOG(INFO) << "Removed duplicated vertices.";

  igl::remove_duplicate_faces_custom(F_);
  LOG(INFO) << "Removed duplicated faces.";

  VectorXd FA;
  igl::doublearea(V_, F_, FA);
  igl::consistent_face_flipping(FA, F_);
  LOG(INFO) << "Flipped faces consistently.";

  LOG(INFO) << "After removing duplicates:";
  LOG(INFO) << " - # vertices: " << n_vertices();
  LOG(INFO) << " - # faces: " << n_faces();

  if (renderer_ != nullptr) {
    renderer_->set_mesh(V_, F_);
  }
}

void LibiglMesh::upsample_mesh(
    const double edge_length_tol, const int max_loop_iters) {
  MatrixXd EL;
  double max_edge_legnth;

  igl::edge_lengths(V_, F_, EL);
  max_edge_legnth = EL.maxCoeff();
  LOG(INFO) << "Before sudivision:";
  LOG(INFO) << " - # vertices: " << n_vertices();
  LOG(INFO) << " - # faces: " << n_faces();
  LOG(INFO) << " - Max edge length: " << max_edge_legnth;

  const int num_iters = std::min(
      (int)std::ceil(std::log2(max_edge_legnth/edge_length_tol)),
      max_loop_iters);
  LOG(INFO) << "# subdivision iterations: " << num_iters;

  LOG(INFO) << "Sudivision started.";
  for (int iter = 0; iter < num_iters; ++iter) {
    LOG(INFO) << " - Iteration: " << iter << " / " << num_iters;
    igl::upsample_non_manifold(V_, F_, edge_length_tol);
  }
  LOG(INFO) << "Sudivision finished.";

  igl::edge_lengths(V_, F_, EL);
  max_edge_legnth = EL.maxCoeff();
  LOG(INFO) << "After sudivision:";
  LOG(INFO) << " - # vertices: " << n_vertices();
  LOG(INFO) << " - # faces: " << n_faces();
  LOG(INFO) << " - Max edge length: " << max_edge_legnth;

  if (renderer_ != nullptr) {
    renderer_->set_mesh(V_, F_);
  }
}

// Last modified: 01-09-2017
void LibiglMesh::processing_subdivide_mesh() {
  // Clean mesh.
  remove_duplicates();

  // Upsample mesh.
  const double kEdgeLengthTol = 0.02;
  const int kMaxLoopIters = 10;
  upsample_mesh(kEdgeLengthTol, kMaxLoopIters);
}
