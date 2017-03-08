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
#include <igl/point_mesh_squared_distance.h>
#include <modules/consistent_face_flippping.h>
#include <modules/remove_small_components.h>
#include <utils/mrf_potts.h>


void LibiglMesh::processing_label_contacting_points(
    const std::string& _out_point_labels_file,
    const double _max_contacting_squared_distance) {
  const int num_points = P_.rows();

  const VectorXi label_set = Utils::unique(FL_);
  MatrixXi P_contacted = MatrixXi::Zero(num_points, label_set.maxCoeff());

  PL_ = VectorXi::Zero(num_points);
  VectorXd min_sqrD = VectorXd::Constant(
      num_points, std::numeric_limits<double>::max());

  // Process for each label mesh.
  for (int i = 0; i < label_set.size(); ++i) {
    const int label = label_set[i];
    // NOTE:
    // Assume that labels start from 1.
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
    igl::point_mesh_squared_distance(P_, label_V, label_F, sqrD, I, C);

    // NOTE:
    // Assume that labels start from 1.
    for (int j = 0; j < num_points; ++j) {
      if (sqrD(j) <= _max_contacting_squared_distance) {
        P_contacted(j, label - 1) = 1;

        if (sqrD(j) < min_sqrD(j)) {
          PL_(j) = label;
          min_sqrD(j) = sqrD(j);
        }
      }
    }
  }

  set_point_label_colors();

  CHECK(Utils::write_eigen_matrix_to_file(_out_point_labels_file, P_contacted));
}
