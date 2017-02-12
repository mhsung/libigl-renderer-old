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


void LibiglMesh::processing_color_barycenter(
    const std::string& _coloring_reference_mesh_file) {
  const double KDefaultMax = 0.3;
  Vector3d ref_bb_min(-KDefaultMax, -KDefaultMax, -KDefaultMax);
  Vector3d ref_bb_max(+KDefaultMax, +KDefaultMax, +KDefaultMax);

  if (_coloring_reference_mesh_file != "") {
    MatrixXd RV;
    MatrixXi RF;
    if (!igl::read_triangle_mesh(_coloring_reference_mesh_file, RV, RF)) {
      LOG(WARNING) << "Can't read the file: '"
                   << _coloring_reference_mesh_file << "'";
      return;
    }
    ref_bb_min = V_.colwise().minCoeff();
    ref_bb_max = V_.colwise().maxCoeff();
  }

  const Vector3d ref_bb_size = (ref_bb_max - ref_bb_min);

  MatrixXd BC;
  igl::barycenter(V_, F_, BC);
  FC_ = MatrixXf::Zero(n_faces(), 3);

  for (int fid = 0; fid < n_faces(); ++fid) {
    RowVector3d center = BC.row(fid);


    RowVector3d color;
    // R[0 ~ 1] <-> y [-KCoordMax, KCoordMax]
    color[0] = (center[1] - ref_bb_min[1]) / ref_bb_size[1];
    // G[0 ~ 1] <-> -x [-KCoordMax, KCoordMax]
    color[1] = (-center[0] - ref_bb_min[0]) / ref_bb_size[0];
    // B[0 ~ 1] <-> abs(z) [0.0, KCoordMax]
    color[2] = std::abs(center[2]) / (0.5 * ref_bb_size[2]);

    for (int i = 0; i < 3; ++i) {
      color[i] = std::min(std::max(color[i], 0.0), 1.0);
    }

    FC_.row(fid) = color.cast<float>();
  }

  if (renderer_ == nullptr) {
    LOG(WARNING) << "Renderer is not set";
  } else {
    renderer_->set_face_colors(FC_);
  }
}
