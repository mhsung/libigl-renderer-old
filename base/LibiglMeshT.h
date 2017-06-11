// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef LIBIGL_MESH_T_H
#define LIBIGL_MESH_T_H

#include <string>
#include <Eigen/Core>
#include <utils/google_tools.h>
#include "LibiglMeshRendererT.h"

using namespace Eigen;

// Declare input variables.
DECLARE_string(mesh);
DECLARE_string(face_labels);
DECLARE_string(point_set);
DECLARE_string(point_labels);
DECLARE_double(azimuth_deg);
DECLARE_double(elevation_deg);
DECLARE_double(theta_deg);
DECLARE_string(projection_matrix);
DECLARE_string(modelview_matrix);
DECLARE_string(snapshot);
DECLARE_string(out_mesh);
DECLARE_string(out_point_set);


class LibiglMeshT {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    LibiglMeshT();
    LibiglMeshT(LibiglMeshRendererT* _renderer);

    int n_vertices() const { return V_.rows(); }
    int n_faces() const { return F_.rows(); }
    int n_points() const { return P_.rows(); }

    bool read_mesh(const std::string& _filename);
    bool read_point_set(const std::string& _filename);

    bool read_face_labels(const std::string& _filename);
    bool write_face_labels(const std::string& _filename);
    void set_face_label_colors();

    bool read_point_labels(const std::string& _filename);
    void set_point_label_colors();

    void update_bounding_box();
    bool write_bounding_box(const std::string& _filename);

    void parse_arguments_and_run();

  protected:
    // Virtual function.
    virtual void processing() = 0;

    void pre_processing();
    void post_processing();

    // Mesh properties.
    std::string mesh_name_;   // Includes extension.

    // Mesh.
    MatrixXd V_;
    MatrixXi F_;
    MatrixXf VC_;
    MatrixXf FC_;
    VectorXi VL_;
    VectorXi FL_;
    MatrixXd VN_;
    MatrixXd FN_;

    // Point set.
    MatrixXd P_;
    MatrixXf PC_;
    VectorXi PL_;
    MatrixXd PN_;

    Vector3d bb_min_;
    Vector3d bb_max_;
    Vector3d center_;
    double radius_;

    // Rendering properties.
    LibiglMeshRendererT* renderer_;
};

#endif	// LIBIGL_MESH_T_H
