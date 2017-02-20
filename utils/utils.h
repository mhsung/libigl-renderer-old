// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <initializer_list>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <Eigen/Core>
#include <Eigen/LU>

#include <utils/google_tools.h>

using namespace Eigen;


namespace Utils {


// -- Random color utils -- //

void random_label_rgb_color(
    const int _label, Vector3f* _color);


// -- Eigen matrix utils -- //

// @_vec: Input vector.
// @_query: Query value.
// @return: Vector indices with query value.
template <typename T>
VectorXi find(const Matrix<T, Dynamic, 1>& _vec, const T _query);

// @_vec: Input vector.
// @return: Unique element vector.
template <typename T>
Matrix<T, Dynamic, 1> unique(const Matrix<T, Dynamic, 1>& _vec);

// @_vec: Input vector.
// @_unique: (Output) Unique element vector.
// @_counts: (Output) Counts of unique elements.
template <typename T>
void count_occurrences(const Matrix<T, Dynamic, 1>& _vec,
                       Matrix<T, Dynamic, 1>* _unique, VectorXi* _counts);

// @_mat: Input matrix.
// @_indices: Row indices to slice.
// @_sub_mat: Sliced sub matrix.
template <typename Scalar, int Rows, int Columns>
Matrix<Scalar, Rows, Columns> slice_rows(
    const Matrix<Scalar, Rows, Columns>& _mat, const VectorXi& _idxs);

// @_mat: Input matrix.
// @_indices: Column indices to slice.
// @_sub_mat: Sliced sub matrix.
template <typename Scalar, int Rows, int Columns>
Matrix<Scalar, Rows, Columns> slice_columns(
    const Matrix<Scalar, Rows, Columns>& _mat, const VectorXi& _idxs);


// -- Eigen I/O utils -- //

template <typename T>
bool write_sequence(
    const std::string& _filepath, const T& _sequence);

template <typename Derived>
bool read_eigen_matrix_from_file(
    const std::string& _filepath,
    MatrixBase<Derived>* _matrix,
    const char _delimiter = ',');

template <typename Derived>
bool write_eigen_matrix_to_file(
    const std::string& _filepath,
    const MatrixBase<Derived>& _matrix,
    const char _delimiter = ',');

template <typename Derived>
bool read_eigen_matrix_from_binary(
    const std::string& _filepath,
    MatrixBase<Derived>* _matrix);

template <typename Derived>
bool write_eigen_matrix_to_binary(
    const std::string& _filepath,
    const MatrixBase<Derived>& _matrix);


// -- Template function implementation -- //

template <typename T>
VectorXi find(const Matrix<T, Dynamic, 1>& _vec, const T _query) {
  VectorXi ret;

  std::list<int> idxs;
  for (int i = 0; i < _vec.size(); ++i) {
    if (_vec[i] == _query) idxs.push_back(i);
  }

  ret = VectorXi(idxs.size());
  int count = 0;
  for (const auto idx : idxs) {
    ret[count++] = idx;
  }

  return ret;
}

template <typename T>
Matrix<T, Dynamic, 1> unique(const Matrix<T, Dynamic, 1>& _vec) {
  std::set<T> elements;
  for (int i = 0; i < _vec.size(); ++i) {
    elements.insert(_vec[i]);
  }

  Matrix<T, Dynamic, 1> ret(elements.size());
  int count = 0;
  for (const auto& element : elements) {
    ret[count++] = element;
  }

  return ret;
}

template <typename T>
void count_occurrences(const Matrix<T, Dynamic, 1>& _vec,
                       Matrix<T, Dynamic, 1>* _unique, VectorXi* _counts) {
  CHECK_NOTNULL(_unique);
  CHECK_NOTNULL(_counts);
  (*_unique) = unique(_vec);

  const int n = _unique->size();
  (*_counts) = VectorXi::Zero(n);

  for (int i = 0; i < _vec.size(); ++i) {
    for (int k = 0; k < n; ++k) {
      if (_vec[i] == (*_unique)[k]) {
        ++((*_counts)[k]);
        break;
      }
    }
  }
}

template <typename Scalar, int Rows, int Columns>
Matrix<Scalar, Rows, Columns> slice_rows(
    const Matrix<Scalar, Rows, Columns>& _mat, const VectorXi& _idxs) {
  const int n = _idxs.size();
  Matrix<Scalar, Rows, Columns> ret;
  ret.resize(n, _mat.cols());

  int count = 0;
  for (int i = 0; i < n; ++i) {
    const int idx = _idxs[i];
    CHECK_LE(idx, _mat.rows());
    ret.row(count) = _mat.row(idx);
    ++count;
  }

  return ret;
}

template <typename Scalar, int Rows, int Columns>
Matrix<Scalar, Rows, Columns> slice_columns(
    const Matrix<Scalar, Rows, Columns>& _mat, const VectorXi& _idxs) {
  const int n = _idxs.size();
  Matrix<Scalar, Rows, Columns> ret;
  ret.resize(_mat.rows(), n);

  int count = 0;
  for (int i = 0; i < n; ++i) {
    const int idx = _idxs[i];
    CHECK_LE(idx, _mat.columns());
    ret.column(count) = _mat.column(idx);
    ++count;
  }

  return ret;
}

template <typename Scalar>
Scalar string_to_number(const std::string& _str) {
  if (_str.size() == 0) return 0;
  std::istringstream sstr(_str);
  Scalar value = 0;
  if (!(sstr >> std::dec >> value)) throw std::invalid_argument(_str);
  return value;
}

template <typename T>
bool write_sequence(
    const std::string& _filepath, const T& _sequence) {
  std::ofstream file(_filepath);
  if (!file.good()) {
    LOG(WARNING) << "Can't write the file: '" << _filepath << "'." << std::endl;
    return false;
  }

  for (auto it = _sequence.begin(); it != _sequence.end(); ++it) {
    file << (*it) << std::endl;
  }

  file.close();
  return true;
}

template <typename Derived>
bool read_eigen_matrix_from_file(
    const std::string& _filepath,
    MatrixBase<Derived>* _matrix,
    const char _delimiter /*= ','*/) {
  CHECK_NOTNULL(_matrix);
  typedef typename Derived::Scalar Scalar;

  std::ifstream file(_filepath);
  if (!file.good()) {
    LOG(WARNING) << "Can't read the file: '" << _filepath << "'";
    return false;
  }

  typedef std::vector<Scalar> StdVector;
  typedef std::unique_ptr<StdVector> StdVectorPtr;
  typedef std::vector<StdVectorPtr> StdMatrix;
  StdMatrix std_matrix;

  std::string line("");
  int num_rows = 0, num_cols = -1;

  for (; std::getline(file, line); ++num_rows) {
    // Stop reading when the line is blank.
    if (line == "") break;
    std::stringstream sstr(line);
    StdVectorPtr vec(new StdVector);

    std::string token("");
    while (std::getline(sstr, token, _delimiter)) {
      // Stop reading when the token is blank.
      if (token == "") break;
      try {
        const Scalar value = string_to_number<Scalar>(token);
        vec->push_back(value);
      }
      catch (std::exception& e) {
        LOG(WARNING) << e.what();
        return false;
      }
    }

    if (num_cols >= 0 && num_cols != vec->size()) {
      LOG(WARNING) << "The number of cols does not match: " <<
                   num_cols << " != " << vec->size();
      return false;
    }

    num_cols = static_cast<int>(vec->size());
    std_matrix.push_back(std::move(vec));
  }

  (*_matrix) = Matrix<Scalar, Dynamic, Dynamic>(
      num_rows, num_cols);
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      (*_matrix)(i, j) = (*std_matrix[i])[j];
    }
  }

  file.close();
  return true;
}

template <typename Derived>
bool write_eigen_matrix_to_file(
    const std::string& _filepath,
    const MatrixBase<Derived>& _matrix,
    const char _delimiter /*= ','*/) {
  std::ofstream file(_filepath);
  if (!file.good()) {
    LOG(WARNING) << "Can't write the file: '" << _filepath << "'";
    return false;
  }

  const IOFormat csv_format(
      FullPrecision, DontAlignCols, std::string(1, _delimiter));
  file << _matrix.format(csv_format);
  file.close();
  return true;
}

template <typename Derived>
bool read_eigen_matrix_from_binary(
    const std::string& _filepath,
    const MatrixBase<Derived>* _matrix) {
  CHECK_NOTNULL(_matrix);

  typedef typename Derived::Scalar Scalar;

  std::ifstream file(_filepath, std::ios::in | std::ios::binary);
  if (!file.good()) {
    LOG(WARNING) << "Can't read the file: '" << _filepath << "'";
    return false;
  }

  int32_t rows = 0, cols = 0;
  file.read((char*) (&rows), sizeof(int32_t));
  file.read((char*) (&cols), sizeof(int32_t));
  _matrix->resize(rows, cols);
  file.read((char*) _matrix->data(), rows * cols * sizeof(Scalar));
  file.close();
  return true;
}

template <typename Derived>
bool write_eigen_matrix_to_binary(
    const std::string& _filepath,
    const MatrixBase<Derived>& _matrix) {
  typedef typename Derived::Scalar Scalar;

  std::ofstream file(_filepath,
                     std::ios::out | std::ios::binary | std::ios::trunc);
  if (!file.good()) {
    LOG(WARNING) << "Can't write the file: '" << _filepath << "'";
    return false;
  }

  int32_t rows = static_cast<int32_t>(_matrix.rows());
  int32_t cols = static_cast<int32_t>(_matrix.cols());
  file.write((char*) (&rows), sizeof(int32_t));
  file.write((char*) (&cols), sizeof(int32_t));
  file.write((char*) _matrix.data(), rows * cols * sizeof(Scalar));
  file.close();
  return true;
}

}

#endif // UTILS_H
