#pragma once

#include "viewer/ViewerTypes.hxx"

namespace viewer {

// Compute the world-space box dimensions for the current view state.
void computeWorldBox(const ViewState &vs, double &worldW, double &worldH);

// Apply the current ViewState to the OpenGL orthographic projection.
void applyOrtho(const ViewState &view);

// Recompute and apply ortho for given framebuffer size.
void resizeAndApplyOrtho(ViewState &view, int fbw, int fbh);

// Pan the view by a screen-pixel delta (right-drag style).
void panView(ViewState &view, double dxPx, double dyPx);

// Zoom the view around a screen-pixel cursor position.
void zoomView(ViewState &view, double scrollDelta, double cursorX, double cursorY);

} // namespace viewer
