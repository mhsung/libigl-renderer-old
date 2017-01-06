// Copyright (C) 2017 Minhyuk Sung <mhsung@cs.stanford.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include "google_tools.h"
#include "LibiglMesh.h"
#ifdef USE_OSMESA
#include "OSMesaMeshRenderer.h"
#else
#include "LibiglMeshViewer.h"
#endif


int main(int argc, char *argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  const int kWidth = 640;
  const int kHeight = 480;

  LibiglMeshRendererT* renderer = nullptr;

#ifdef USE_OSMESA
  OSMesaMeshRenderer osmesa_renderer(kWidth, kHeight);
  renderer = &osmesa_renderer;
#else
  LibiglMeshViewer libigl_renderer(kWidth, kHeight);
  renderer = &libigl_renderer;
#endif

  LibiglMesh mesh(renderer);
  mesh.parse_arguments_and_run();
}
