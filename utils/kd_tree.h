// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef KD_TREE_H
#define KD_TREE_H

#include <memory>
#include <Eigen/Core>
#include <utils/google_tools.h>
#include <utils/nanoflann/nanoflann.hpp>

using namespace Eigen;


namespace KdTree {

/*
 * Example:
 * // D = Data points, Q = Query points (in each row).
 * // NOTE: Do not delete D before retrieving closest points.
 * MatrixXd D, Q;
 * KdTree::::KdTreeXd tree(P);
 * const VectorXi ret = KdTree::::find_closest_points(tree, Q);
 */

// Reference:
// https://github.com/jlblancoc/nanoflann/blob/master/tests/test_main.cpp

template<typename Scalar,int Rows,int Cols>
struct KdTree {
  // @_data_pts: N-dimensional data points (each row, N = 2, 3, or 4).
  KdTree(const Matrix<Scalar, Dynamic, Cols>& _data_pts)
      : data_pts_(_data_pts),
        adaptor_(_data_pts.cols(), *this,
                 nanoflann::KDTreeSingleIndexAdaptorParams(10)) {
    CHECK(dim() == 2 || dim() == 3 || dim() == 4)
    << "Dimension must be 2, 3, or 4.";
    adaptor_.buildIndex();
  }

  int dim() const { return data_pts_.cols(); }

  typedef nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<Scalar, KdTree<Scalar, Rows, Cols> >,
      KdTree<Scalar, Rows, Cols>> KdTreeAdaptor;

  const Matrix<Scalar, Rows, Cols>& data_pts_;
  KdTreeAdaptor adaptor_;

  // Template functions required in nanoflann.
  inline size_t kdtree_get_point_count() const { return data_pts_.rows(); }

  inline Scalar kdtree_distance(
      const Scalar* p1, const size_t idx2, size_t /* size*/) const {
    const auto& p2 = data_pts_.row(idx2);
    Scalar ret = 0;
    for (int i = 0; i < dim(); ++i) {
      const Scalar diff = p1[i] - p2[i];
      ret += (diff * diff);
    }
    return ret;
  }

  inline Scalar kdtree_get_pt(const size_t idx, int dim) const {
    const auto& p = data_pts_.row(idx);
    return p[dim];
  }

  template<class BBOX>
  bool kdtree_get_bbox(BBOX& /* bb*/ ) const { return false; }
};

typedef KdTree<float, Dynamic, Dynamic> KdTreeXf;
typedef KdTree<double, Dynamic, Dynamic> KdTreeXd;


// @_tree: k-d Tree.
// @_query_pts: N-dimensional query points (each row, N = 2, 3, or 4).
// @return: Closest point indices.
template<typename Scalar,int Rows1,int Cols1,int Rows2,int Cols2>
VectorXi find_closest_points(
    KdTree<Scalar, Rows1, Cols1>& _tree,
    const Matrix<Scalar, Rows2, Cols2>& _query_pts) {
  CHECK_EQ(_query_pts.cols(), _tree.dim())
  << "Data and query point dimensions do not match.";

  const int num_query = _query_pts.rows();
  VectorXi ret(num_query);

  for (int i = 0; i < num_query; ++i) {
    // Find the closest data point.
    const Matrix<Scalar, 1, Dynamic> query_pt = _query_pts.row(i);
    size_t ret_idx;
    Scalar sqr_dist;
    nanoflann::KNNResultSet<Scalar> results(1);
    results.init(&ret_idx, &sqr_dist);
    CHECK(_tree.adaptor_.findNeighbors(
        results, query_pt.data(), nanoflann::SearchParams(10)));
    ret[i] = ret_idx;
  }

  return ret;
}

}

#endif // KD_TREE_H
