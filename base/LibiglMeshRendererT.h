// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef LIBIGL_MESH_RENDERER_T_H
#define LIBIGL_MESH_RENDERER_T_H

#include <string>
#include <vector>
#include <Eigen/Core>

using namespace Eigen;


class LibiglMeshRendererT {
  public:
    LibiglMeshRendererT(const int _width, const int _height);
    virtual ~LibiglMeshRendererT() {};

    bool read_projection(const std::string& _filename);
    bool read_modelview(const std::string& _filename);
    bool write_projection(const std::string& _filename) const;
    bool write_modelview(const std::string& _filename) const;
    void set_scene_pos(const Vector3f& _center, const float _radius);

    // Camera parameters: azimuch, elevation, theta (in degrees).
    Vector3f get_camera_params() const;
    void set_camera_params(
        const Vector3f& camera_params, const float _radius = 1.0f);

    // Virtual functions.
    virtual const Matrix4f& get_projection() const = 0;
    virtual const Matrix4f& get_modelview() const = 0;
    virtual void set_projection(const Matrix4f& _projection) = 0;
    virtual void set_modelview(const Matrix4f& _modelview) = 0;

    virtual void set_mesh(
        const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F) = 0;
    virtual void set_vertex_colors(const Eigen::MatrixXf& _VC) = 0;
    virtual void set_face_colors(const Eigen::MatrixXf& _FC) = 0;

    virtual void set_points(const Eigen::MatrixXd& _P) = 0;
    virtual void set_point_colors(const Eigen::MatrixXf& _PC) = 0;

    virtual void set_keypoints(const Eigen::MatrixXd& _key_P) = 0;
    virtual void set_keypoint_colors(const Eigen::MatrixXf& _key_PC) = 0;

    virtual void run_loop() = 0;
    virtual bool snapshot(const std::string& _filename) = 0;

  protected:
    const int width_;
    const int height_;
};

#endif	// LIBIGL_MESH_RENDERER_T_H

