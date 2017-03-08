// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef LIBIGL_MESH_H
#define LIBIGL_MESH_H

#include "LibiglMeshT.h"

using namespace Eigen;


class LibiglMesh : public LibiglMeshT {
  public:
    LibiglMesh() : LibiglMeshT() {};
    LibiglMesh(LibiglMeshRendererT* _renderer) : LibiglMeshT(_renderer) {};

    void remove_duplicates();

    void upsample_mesh(
        const double edge_length_tol, const int max_loop_iters = 10);

    void transform_mesh(const std::string& _filename);

    // Move center to (0,0,0) and scale to bounding box diagonal 1.
    void normalize_mesh(MatrixXd& V);


  protected:
    virtual void processing();

    void processing_subdivide_mesh();

    void processing_project_pts_labels_to_mesh(
        const std::string& _point_set_file,
        const std::string& _point_labels_file,
        const std::string& _out_face_labels_file);

    void processing_disassemble_to_parts(
        const std::string& _out_part_mesh_dir,
        const std::string& _out_part_mesh_unnormalized_dir,
        const std::string& _out_part_mesh_face_map_dir);

    void processing_disassemble_to_components(
        const std::string& _out_component_mesh_dir,
        const std::string& _out_component_mesh_unnormalized_dir,
        const std::string& _out_component_mesh_face_map_file,
        const int _min_num_components,
        const int _min_component_bbox_diagonal,
        const bool _find_symmetry = false);

    void processing_sample_points(
        const int _num_sample_points,
        const std::string& _out_point_set_dir,
        const std::string& _out_pca_transformation_dir,
        const std::string& _out_position_dir,
        const std::string& _out_centered_point_set_dir,
        const std::string& _out_pca_aligned_point_set_dir);

    void processing_color_barycenter(
        const std::string& _coloring_reference_mesh_file);

    void processing_label_contacting_points(
        const std::string& _point_set_file,
        const std::string& _out_point_labels_file,
        const double _max_contacting_squared_distance);
};

#endif	// LIBIGL_MESH_H
