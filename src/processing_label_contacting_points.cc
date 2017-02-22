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


void LibiglMesh::processing_label_contacting_points(
    const std::string& _point_set_file,
    const std::string& _out_point_labels_file,
    const double _max_contacting_squared_distance) {
  // Read point set and labels.
  MatrixXd P;
  if (!read_point_set(_point_set_file, &P)) return;
  const int num_points = P.rows();

  const VectorXi label_set = Utils::unique(FL_);
  MatrixXi PL = MatrixXi::Zero(num_points, label_set.maxCoeff() + 1);

  // Process for each label mesh.
  for (int i = 0; i < label_set.size(); ++i) {
    const int label = label_set[i];
    // NOTE:
    // Assume that there is no zero label.
    CHECK_GT(label, 0);

    const VectorXi label_fids = Utils::find(FL_, label);
    const MatrixXi label_F_old = Utils::slice_rows(F_, label_fids);
    MatrixXd label_V;
    MatrixXi label_F;
    VectorXi IX;
    igl::remove_unreferenced(V_, label_F_old, label_V, label_F, IX);

    VectorXd sqrD;
    VectorXi I;
    MatrixXd C;
    igl::point_mesh_squared_distance(P, label_V, label_F, sqrD, I, C);

    for (int j = 0; j < num_points; ++j) {
      if (sqrD(j) <= _max_contacting_squared_distance) PL(j, label) = 1;
    }
  }

  LOG(INFO) << label_set.transpose();
  LOG(INFO) << PL.colwise().sum();

  CHECK(Utils::write_eigen_matrix_to_file(_out_point_labels_file, PL));
}
