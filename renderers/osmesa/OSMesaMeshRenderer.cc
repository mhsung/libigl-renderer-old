// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "OSMesaMeshRenderer.h"

#include <iostream>
#include <igl/material_colors.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/png/writePNG.h>
#include <utils/google_tools.h>
#include <utils/glut_geometry.h>


OSMesaMeshRenderer::OSMesaMeshRenderer(const int _width, const int _height)
  : LibiglMeshRendererT(_width, _height) {
  projection_ = Matrix4f::Identity();
  modelview_ = Matrix4f::Identity();

  initialize_osmesa();
  initialize_opengl();
}

OSMesaMeshRenderer::~OSMesaMeshRenderer() {
  free(frame_buffer_);
  OSMesaDestroyContext(context_);
}

const Matrix4f& OSMesaMeshRenderer::get_projection() const {
  return projection_;
}

const Matrix4f& OSMesaMeshRenderer::get_modelview() const {
  return modelview_;
}

void OSMesaMeshRenderer::set_projection(const Matrix4f& _projection) {
  projection_ = _projection;
}

void OSMesaMeshRenderer::set_modelview(const Matrix4f& _modelview) {
  modelview_ = _modelview;
}

void OSMesaMeshRenderer::set_mesh(
    const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F) {
  CHECK_EQ(_V.cols(), 3);
  V_.resize(_V.rows(), 3);
  for (int i = 0 ; i < _V.rows(); ++i) V_.row(i) = _V.row(i);
  F_.resize(_F.rows(), 3);
  for (int i = 0 ; i < _F.rows(); ++i) F_.row(i) = _F.row(i);

  igl::per_face_normals(V_, F_, FN_);
  igl::per_vertex_normals(V_, F_, FN_, VN_);
}

void OSMesaMeshRenderer::set_vertex_colors(const Eigen::MatrixXf& _VC) {
  CHECK_EQ(_VC.cols(), 3);
  VC_.resize(_VC.rows(), 3);
  for (int i = 0 ; i < _VC.rows(); ++i) VC_.row(i) = _VC.row(i);
}

void OSMesaMeshRenderer::set_face_colors(const Eigen::MatrixXf& _FC) {
  CHECK_EQ(_FC.cols(), 3);
  FC_.resize(_FC.rows(), 3);
  for (int i = 0 ; i < _FC.rows(); ++i) FC_.row(i) = _FC.row(i);
}

void OSMesaMeshRenderer::set_points(const Eigen::MatrixXd& _P) {
  CHECK_EQ(_P.cols(), 3);
  P_.resize(_P.rows(), 3);
  for (int i = 0 ; i < _P.rows(); ++i) P_.row(i) = _P.row(i);
}

void OSMesaMeshRenderer::set_point_colors(const Eigen::MatrixXf& _PC) {
  CHECK_EQ(_PC.cols(), 3);
  PC_.resize(_PC.rows(), 3);
  for (int i = 0 ; i < _PC.rows(); ++i) PC_.row(i) = _PC.row(i);
}

void OSMesaMeshRenderer::set_keypoints(const Eigen::MatrixXd& _key_P) {
  CHECK_EQ(_key_P.cols(), 3);
  key_P_.resize(_key_P.rows(), 3);
  for (int i = 0 ; i < _key_P.rows(); ++i) key_P_.row(i) = _key_P.row(i);
}

void OSMesaMeshRenderer::set_keypoint_colors(const Eigen::MatrixXf& _key_PC) {
  CHECK_EQ(_key_PC.cols(), 3);
  key_PC_.resize(_key_PC.rows(), 3);
  for (int i = 0 ; i < _key_PC.rows(); ++i) key_PC_.row(i) = _key_PC.row(i);
}

void OSMesaMeshRenderer::run_loop() {
}

bool OSMesaMeshRenderer::snapshot(const std::string& _filename) {
  // Draw mesh.
  render();

  // Allocate temporary buffers.
  Matrix<unsigned char, Dynamic, Dynamic> R(width_, height_);
  Matrix<unsigned char, Dynamic, Dynamic> G(width_, height_);
  Matrix<unsigned char, Dynamic, Dynamic> B(width_, height_);
  Matrix<unsigned char, Dynamic, Dynamic> A(width_, height_);

  const unsigned char *ptr = static_cast<unsigned char *>(frame_buffer_);
  for (int y = height_ - 1; y >= 0; --y) {
    for (int x = 0; x < width_; ++x) {
      const int i = (y * width_ + x) * 4;
      R(x, y) = ptr[i];
      G(x, y) = ptr[i + 1];
      B(x, y) = ptr[i + 2];
      A(x, y) = ptr[i + 3];
    }
  }

  // Save it to a PNG.
  return igl::png::writePNG(R, G, B, A, _filename + std::string(".png"));
}

void OSMesaMeshRenderer::initialize_osmesa() {
  context_ = OSMesaCreateContextExt(OSMESA_RGBA, 16, 0, 0, NULL);
  if (!context_) {
    LOG(ERROR) << "OSMesaCreateContext failed!";
  }

  frame_buffer_ = malloc(width_ * height_ * 4 * sizeof(GLubyte));
  if (!frame_buffer_) {
    LOG(ERROR) << "Allocating image buffer failed!";
  }

  if (!OSMesaMakeCurrent(context_, frame_buffer_, GL_UNSIGNED_BYTE,
        width_, height_)) {
    LOG(ERROR) << "OSMesaMakeCurrent failed!";
  }

  int z, s, a;
  glGetIntegerv(GL_DEPTH_BITS, &z);
  glGetIntegerv(GL_STENCIL_BITS, &s);
  glGetIntegerv(GL_ACCUM_RED_BITS, &a);
  LOG(INFO) << "Depth = " << z << "  Stencil = " << s << "  Accum" << a;
}

void OSMesaMeshRenderer::initialize_opengl() {
  // OpenGL state
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glDisable(GL_DITHER);
  glEnable(GL_DEPTH_TEST);

  // Material
  set_default_material();

  // Lighting
  set_default_light();

  // Fog
  GLfloat fogColor[4] = { 0.3, 0.3, 0.4, 1.0 };
  glFogi(GL_FOG_MODE, GL_LINEAR);
  glFogfv(GL_FOG_COLOR, fogColor);
  glFogf(GL_FOG_DENSITY, 0.35);
  glHint(GL_FOG_HINT, GL_DONT_CARE);
  glFogf(GL_FOG_START, 5.0f);
  glFogf(GL_FOG_END, 25.0f);

  // Scene pos and size
  set_scene_pos(Vector3f::Zero(), 1.0f);
}


void OSMesaMeshRenderer::set_default_material() {
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, igl::SILVER_AMBIENT);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, igl::SILVER_DIFFUSE);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, igl::SILVER_SPECULAR);
}

void OSMesaMeshRenderer::set_default_light() {
  glLoadIdentity();

  // Front light
  GLfloat pos1[] = { 1, 1, -0.2, 0.0 };
  GLfloat pos2[] = { -1, 1, -0.2, 0.0 };
  GLfloat pos3[] = { 0, 0, 1, 0.0 };
  //GLfloat col1[] = { 0.7, 0.7, 0.8, 1.0 };
  //GLfloat col2[] = { 0.8, 0.7, 0.7, 1.0 };
  //GLfloat col3[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat col1[] = { 0.21, 0.21, 0.24, 1.0 };
  GLfloat col2[] = { 0.24, 0.21, 0.21, 1.0 };
  GLfloat col3[] = { 0.3, 0.3, 0.3, 1.0 };

  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_POSITION, pos1);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, col1);
  glLightfv(GL_LIGHT0, GL_SPECULAR, col1);

  glEnable(GL_LIGHT1);
  glLightfv(GL_LIGHT1, GL_POSITION, pos2);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, col2);
  glLightfv(GL_LIGHT1, GL_SPECULAR, col2);

  glEnable(GL_LIGHT2);
  glLightfv(GL_LIGHT2, GL_POSITION, pos3);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, col3);
  glLightfv(GL_LIGHT2, GL_SPECULAR, col3);
}

void OSMesaMeshRenderer::render() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadMatrixf(projection_.data());
  glMatrixMode(GL_MODELVIEW);
  glLoadMatrixf(modelview_.data());

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  //glDisable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);

  glCullFace(GL_FRONT_AND_BACK);
  glEnable(GL_NORMALIZE);

  render_point_set();
  render_keypoints();
  render_mesh();

  set_default_material();
}

void OSMesaMeshRenderer::render_mesh() {
  const bool has_vertex_color = (VC_.rows() == V_.rows());
  const bool has_face_color = (FC_.rows() == F_.rows());

  if (has_vertex_color && has_face_color) {
    LOG(ERROR) << "Both vertices and faces have colors.";
  }

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_DOUBLE, 0, V_.data());

  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_DOUBLE, 0, VN_.data());

  glBegin(GL_TRIANGLES);
  for (int fid = 0; fid < F_.rows(); ++fid) {
    if (has_face_color) {
      const Vector4f color(FC_(fid, 0), FC_(fid, 1), FC_(fid, 2), 1.0f);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color.data());
    }

    for (int i = 0; i < 3; ++i) {
      const int vid = F_(fid, i);
      if (has_vertex_color) {
        const Vector4f color(VC_(vid, 0), VC_(vid, 1), VC_(vid, 2), 1.0f);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color.data());
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color.data());
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color.data());
      }
      glArrayElement(vid);
    }
  }
  glEnd();

  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
}

/*
void OSMesaMeshRenderer::render_point_set() {
  const bool has_point_color = (PC_.rows() == P_.rows());

  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_DOUBLE, 0, P_.data());

  glPointSize(4.0);
  glBegin(GL_POINTS);
  for (int pid = 0; pid < P_.rows(); ++pid) {
    if (has_point_color) {
      const Vector3f color = PC_.row(pid);
      glColor3fv(color.data());
    } else {
      glColor4fv(igl::MAYA_GREY.data());
    }
    glArrayElement(pid);
  }
  glEnd();

  glDisableClientState(GL_VERTEX_ARRAY);
}
*/

void OSMesaMeshRenderer::render_point_set() {
  const bool has_point_color = (PC_.rows() == P_.rows());
  const double radius = 0.01;

  for (int pid = 0; pid < P_.rows(); ++pid) {
    if (has_point_color) {
      //const Vector3f color = PC_.row(pid);
      //glColor3fv(color.data());
      Vector4f color;
      for (int i = 0; i < 3; ++i) color[i] = PC_.row(pid)[i];
      color[3] = 1.0f;
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color.data());
    } else {
      set_default_material();
    }

    glPushMatrix();
    glTranslated(P_(pid, 0), P_(pid, 1), P_(pid, 2));
    glutSolidSphere(radius, 16, 16);
    glPopMatrix();
  }

  set_default_material();
}

void OSMesaMeshRenderer::render_keypoints() {
  const bool has_keypoint_color = (key_PC_.rows() == key_P_.rows());
  const double radius = 0.04;

  for (int pid = 0; pid < key_P_.rows(); ++pid) {
    if (has_keypoint_color) {
      //const Vector3f color = key_PC_.row(pid);
      //glColor3fv(color.data());
      Vector4f color;
      for (int i = 0; i < 3; ++i) color[i] = key_PC_.row(pid)[i];
      color[3] = 1.0f;
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color.data());
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color.data());
    } else {
      set_default_material();
    }

    glPushMatrix();
    glTranslated(key_P_(pid, 0), key_P_(pid, 1), key_P_(pid, 2));
    glutSolidSphere(radius, 16, 16);
    glPopMatrix();
  }

  set_default_material();
}

