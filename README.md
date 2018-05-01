libigl-renderer
=================
Author: Minhyuk Sung<br>
libigl-based renderer.<br>
libigl: http://libigl.github.io/libigl/<br>
OSMesa rendering is supported. Set '-DUSE_OSMESA=TRUE' in cmake to use OSMesa renderer.<br>
Copy 'include/LibiglMesh.h.blank' to 'include/LibiglMesh.h' and add functions.

### Requirements
sudo apt-get install libosmesa6-dev libgflags-dev libgoogle-glog-dev libhdf5-serial-dev

### Build
mkdir build && cd build
cmake .. -DUSE_OSMESA=TRUE
make -j4

