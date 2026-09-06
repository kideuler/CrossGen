#pragma once

// Export.hxx -- the figure side of the viewer: a light background, and a way
// to get what is on screen out as a file a paper can include.
//
// Two things here, and they are coupled. The viewer's colours were all chosen
// against a 0.1 grey; a paper wants white, and a figure has to be the thing
// that was on screen rather than a second drawing of it. So the background and
// the colours move together: setLightBackground() flips the clear colour and
// puts every colour the render code asks for through adaptForLight(), which is
// why the draw calls go through color3f()/color4f() here instead of glColor3f
// directly. Nothing in Render.cxx has to know which background it is drawing
// against.

#include <array>
#include <cstddef>
#include <string>

namespace viewer {

// ── theme ────────────────────────────────────────────────────────────────────

// White background (the default, since that is what figures want) or the dark
// one the viewer was written against.
void setLightBackground(bool on);
bool lightBackground();

// The clear colour for the current theme.
std::array<float, 3> backgroundColor();

// Remap a colour chosen for a dark background so it reads on a light one; a
// no-op when the dark theme is active.
//
// Two cases, because two things go wrong. A near-neutral has its lightness
// inverted -- the white an interface edge is drawn in becomes near-black, a
// mid grey stays a mid grey, since that already reads either way. A saturated
// hue keeps its hue: it is pulled away from white by dropping its darkest
// channel (on white a colour reads by how much it is *not* white, so a pale
// hue reads as nothing at all) and then capped in luminance, which is what
// rescues the greens and yellows that are bright by nature.
void adaptForLight(float &r, float &g, float &b);

// Drop-in replacements for glColor3f/glColor4f that apply the theme. Every
// colour the viewer sets goes through these. The ...Raw forms skip the
// adaptation, for a colour that is already right against either background.
void color3f(float r, float g, float b);
void color4f(float r, float g, float b, float a);
void color3fRaw(float r, float g, float b);
void color4fRaw(float r, float g, float b, float a);

// ── render scale ─────────────────────────────────────────────────────────────

// Line widths, point sizes and the bitmap font are all in pixels, so a
// supersampled raster export would draw them at a fraction of the weight they
// have on screen -- a 4x image with 1px edges is a picture of a different
// mesh. Every such call goes through the wrappers below and is multiplied by
// this scale, which the PNG export sets to its supersampling factor.
void  setRenderScale(float s);
float renderScale();

// glLineWidth/glPointSize with the render scale applied, and -- while a vector
// capture is running -- with a marker pushed into the feedback stream, since
// the stream itself carries no width. Clamped to what the driver will accept.
void lineWidth(float w);
void pointSize(float s);

// The scale the bitmap font should draw at. Render.cxx multiplies its own
// pixel scale by this.
float textScale();

// ── chrome ───────────────────────────────────────────────────────────────────

// The coordinate axis and its tick labels. It is drawn from fifteen places, so
// the switch lives inside drawAxis() rather than at every call site. Off with
// the rest of the chrome for a figure: a green crosshair ruled across the model
// is a thing the viewer wants and a paper does not.
void setAxisVisible(bool on);
bool axisVisible();

// ── vector capture (SVG) ─────────────────────────────────────────────────────

// One frame captured as SVG through the GL feedback buffer: the scene is drawn
// a second time with rasterisation switched off, and what comes back is the
// primitive stream in draw order, already through the projection and the
// viewport, so it is in window coordinates. Depth testing is off in this
// viewer and everything is painted back to front, so stream order is paint
// order and there is nothing to sort -- which is the whole reason this is
// ~300 lines rather than a dependency on gl2ps.
//
// Usage, with the context current:
//
//     beginVectorCapture(capacity);
//     drawScene();                       // exactly as for the screen
//     writeVectorCapture(path, w, h);
//
// The capacity is in floats and cannot be grown once the capture has started
// (the primitives are already gone), so an Overflow return means: draw it
// again with a bigger one. exportSvg() in CrossGenWidget does that loop.
enum class CaptureResult {
    Ok,
    Overflow,    // the feedback buffer was too small; retry with more
    Empty,       // nothing was drawn, or no capture was running
    WriteFailed, // the file could not be opened
};

bool beginVectorCapture(std::size_t capacityFloats);

// The captured frame as an SVG document in memory. This is the form the PDF
// export wants: the figure is written by handing the document to a renderer
// that paints it onto a vector page, so nothing has to go through a file on
// the way. Ends the capture, exactly as writeVectorCapture does.
CaptureResult buildVectorCapture(std::string &svg, int fbw, int fbh);

// The same document, straight to a file.
CaptureResult writeVectorCapture(const std::string &path, int fbw, int fbh);
bool capturing();

// How many primitives the last capture turned into SVG elements, for the
// message the viewer logs.
struct CaptureStats {
    std::size_t polygons{0};
    std::size_t lines{0};
    std::size_t points{0};
    std::size_t elements{0}; // after batching
    std::size_t floatsUsed{0};
};
const CaptureStats &lastCaptureStats();

} // namespace viewer
