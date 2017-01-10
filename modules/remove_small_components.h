// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_REMOVE_SMALL_COMPONENTS_H
#define IGL_REMOVE_SMALL_COMPONENTS_H

#include <igl/igl_inline.h>
#include <igl/doublearea.h>
#include <igl/facet_components.h>
#include <igl/unique.h>

#include <map>
#include <vector>
#include <Eigen/Core>

namespace igl
{
  template <
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE void remove_small_components(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedF>& newF,
    const float eps = 1.0E-3f);

  template <
    typename MatV,
    typename MatF>
  IGL_INLINE void remove_small_components(
    MatV& V,
    MatF& F,
    const float eps = 1.0E-3f);
}

//#ifndef IGL_STATIC_LIBRARY
//#  include "remove_small_components.cpp"
//#endif


template <
    typename DerivedV,
    typename DerivedF>
IGL_INLINE void igl::remove_small_components(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedF>& newF,
    const float eps)
{
  typedef typename DerivedV::Scalar ScalarV;

  // Find components.
  Eigen::VectorXi C;
  igl::facet_components(F, C);

  Eigen::VectorXi C_set;
  igl::unique(C, C_set);
  std::map<int, int> C_to_idx;
  const int num_C_set = C_set.size();
  for (int cid = 0; cid < num_C_set; ++cid) {
    C_to_idx[ C_set[cid] ] = cid;
  }

  // Compute face areas.
  Eigen::Matrix<ScalarV, Eigen::Dynamic, 1> DblFA;
  igl::doublearea(V, F, DblFA);

  // Compute component areas.
  Eigen::Matrix<ScalarV, Eigen::Dynamic, 1> DblCA(num_C_set);
  DblCA.setZero();
  const int nF = F.rows();
  for (int fid = 0; fid < nF; ++fid) {
    const int cid = C_to_idx[ C[fid] ];
    DblCA[cid] += DblFA[fid];
  }

  // Remove faces in small components.
  std::vector<bool> is_component_small(num_C_set, false);
  for (int cid = 0; cid < num_C_set; ++cid) {
    if (DblCA[cid] < 2 * eps) is_component_small[cid] = true;
  }

  std::vector<int> new_to_old_fids;
  new_to_old_fids.reserve(nF);
  for (int fid = 0; fid < nF; ++fid) {
    const int cid = C_to_idx[ C[fid] ];
    if (!is_component_small[cid]) new_to_old_fids.push_back(fid);
  }

  const int new_nF = new_to_old_fids.size();
  if (new_nF == nF) {
    newF = F;
    return;
  }

  newF.resize(new_nF, F.cols());
  for (int fid = 0; fid < new_nF; ++fid) {
    newF.row(fid) = F.row(new_to_old_fids[fid]);
  }
}

template <
    typename MatV,
    typename MatF>
IGL_INLINE void igl::remove_small_components(
    MatV& V,
    MatF& F,
    const float eps)
{
  const MatF F_copy = F;
  remove_small_components(V,F_copy,F,eps);
}

#endif
