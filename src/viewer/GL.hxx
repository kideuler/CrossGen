#pragma once

// OpenGL includes for the viewer (Qt6 / platform-native).
// GLFW has been replaced by Qt6; OpenGL context is managed by QOpenGLWidget.

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif
