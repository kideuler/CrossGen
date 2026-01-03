#pragma once

// OpenGL/GLFW includes for the viewer.
// We centralize them to avoid including GL headers before GLFW when using GLFW_INCLUDE_NONE.

#include <GLFW/glfw3.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
