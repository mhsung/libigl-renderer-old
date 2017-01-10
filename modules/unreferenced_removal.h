// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_UNREFERENCED_REMOVAL_H
#define IGL_UNREFERENCED_REMOVAL_H

#include <igl/igl_inline.h>
#include <igl/sort.h>
#include <utils/kd_tree.h>

#include <Eigen/Core>

namespace igl
{
  template <
    typename DerivedV,
    typename DerivedF>
  IGL_INLINE void unreferenced_vertex_removal(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedV>& newV,
    Eigen::PlainObjectBase<DerivedF>& newF);

  template <
    typename MatV,
    typename MatF>
  IGL_INLINE void unreferenced_vertex_removal(
    MatV& V,
    MatF& F);
}

//#ifndef IGL_STATIC_LIBRARY
//#  include "unreferenced_removal.cpp"
//#endif


template <
    typename DerivedV,
    typename DerivedF>
IGL_INLINE void igl::unreferenced_vertex_removal(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedV>& newV,
    Eigen::PlainObjectBase<DerivedF>& newF)
{
  const int nV = V.rows();
  const int nF = F.rows();

  // Find unreferenced vertices.
  std::vector<bool> is_vertex_referenced(nV, false);
  int num_vertex_referenced = 0;
  for (int fid = 0; fid < nF; ++fid) {
    for (int i = 0; i < F.cols(); ++i) {
      if (!is_vertex_referenced[F(fid, i)]) {
        ++num_vertex_referenced;
        is_vertex_referenced[F(fid, i)] = true;
      }
    }
  }

  newV.resize(nV, V.cols());
  std::vector<int> old_to_new_vids(nV, -1);
  int n_vid = 0;
  for (int vid = 0; vid < nV; ++vid) {
    if (is_vertex_referenced[vid]) {
      newV.row(n_vid) = V.row(vid);
      old_to_new_vids[vid] = n_vid;
      ++n_vid;
    }
  }

  newF.resize(F.rows(), F.cols());
  for (int fid = 0; fid < F.rows(); ++fid) {
    newF.row(fid) << old_to_new_vids[ F.row(fid)[0] ],
        old_to_new_vids[ F.row(fid)[1] ],
        old_to_new_vids[ F.row(fid)[2] ];
  }
}

template <
    typename MatV,
    typename MatF>
IGL_INLINE void igl::unreferenced_vertex_removal(
    MatV& V,
    MatF& F)
{
  const MatV V_copy = V;
  const MatF F_copy = F;
  unreferenced_vertex_removal(V_copy,F_copy,V,F);
}

#endif
