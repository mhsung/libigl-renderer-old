// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PCA_H
#define IGL_PCA_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>

namespace igl
{
  // Run PCA.
  // P: n x k points.
  // T: Affine transformation matrix.
  // Aligned_P_transpose = T * P_transpose.
  //
  template <typename Scalar, int RowP, int ColumnP>
  IGL_INLINE void PCA(
    const Eigen::Matrix<Scalar, RowP, ColumnP>& P,
    Eigen::Transform<Scalar, 3, Eigen::Affine>& T);
}

//#ifndef IGL_STATIC_LIBRARY
//#  include "PCA.cpp"
//#endif


template <typename Scalar, int RowP, int ColumnP>
IGL_INLINE void igl::PCA(
    const Eigen::Matrix<Scalar, RowP, ColumnP>& P,
    Eigen::Transform<Scalar, 3, Eigen::Affine>& T)
{
  assert(P.cols() == 3);
  using namespace Eigen;

  const Matrix<Scalar, ColumnP, RowP> P_t = P.transpose();
  const Matrix<Scalar, 3, 1>  t_t(P_t.rowwise().mean());

  const Matrix<Scalar, ColumnP, RowP> P_t_centered = P_t.colwise() - t_t;
  JacobiSVD<Matrix<Scalar, ColumnP, RowP>> svd(P_t_centered, ComputeFullU);

  // NOTE:
  // Singular values are always sorted in decreasing order.
  // https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html
  Matrix<Scalar, 3, 3> R_t = svd.matrixU();
  if ((R_t.col(0).cross(R_t.col(1))).dot(R_t.col(2)) < 0)
    R_t.col(2) = -R_t.col(2);

  T.pretranslate(-t_t);
  T.prerotate(R_t.inverse());
};

#endif

