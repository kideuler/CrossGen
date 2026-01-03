#pragma once

#include "viewer/ViewerTypes.hxx"

#include "viewer/GL.hxx"

namespace viewer {

// Setup callbacks for pan/zoom and framebuffer resize.
// The caller must set the GLFW window user pointer to a `ViewState*` before calling.
void installInteractionCallbacks(GLFWwindow *window);

// Apply the current ViewState to the OpenGL orthographic projection.
void applyOrtho(const ViewState &view);

// Recompute and apply ortho for given framebuffer size.
void resizeAndApplyOrtho(ViewState &view, int fbw, int fbh);

} // namespace viewer
