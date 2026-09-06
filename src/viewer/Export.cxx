#include "viewer/Export.hxx"
#include "viewer/GL.hxx"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace viewer {
namespace {

bool  g_light = true;   // white by default: figures are what this is for
float g_scale = 1.0f;
bool  g_axis  = true;

// ── feedback capture state ───────────────────────────────────────────────────
bool                 g_capturing = false;
std::vector<GLfloat> g_feedback;
CaptureStats         g_stats;

// Markers pushed into the feedback stream by lineWidth()/pointSize(). A
// passthrough token carries exactly one float, so the value has to say which
// of the two it is; the offsets are far outside any plausible pixel width, and
// the width itself is the fractional remainder.
constexpr float kLineWidthMark = 100000.0f;
constexpr float kPointSizeMark = 200000.0f;

// GL_3D_COLOR in an RGBA context: x, y, z, then r, g, b, a.
constexpr int kVertexFloats = 7;

// What the driver will actually accept for glLineWidth. Asking for 14 px on a
// 4x supersampled export and silently getting 10 is fine on screen but would
// be a lie in the SVG, so the clamp applies to the GL call only -- the capture
// records the width that was asked for.
float maxGLLineWidth() {
    static float cached = -1.0f;
    if (cached < 0.0f) {
        GLfloat range[2] = {1.0f, 1.0f};
        glGetFloatv(GL_LINE_WIDTH_RANGE, range);
        cached = (range[1] >= 1.0f) ? range[1] : 1.0f;
    }
    return cached;
}

// ── SVG helpers ──────────────────────────────────────────────────────────────

// Two decimals is well under a pixel and keeps the file a fraction of the size
// printf's default six would make it.
void appendNum(std::string &out, float v) {
    if (!std::isfinite(v)) v = 0.0f;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.2f", static_cast<double>(v));
    // Trim the trailing zeros and a bare '.' that %.2f leaves on round numbers.
    std::size_t n = std::strlen(buf);
    if (std::strchr(buf, '.') != nullptr) {
        while (n > 0 && buf[n - 1] == '0') --n;
        if (n > 0 && buf[n - 1] == '.') --n;
    }
    out.append(buf, n);
}

void appendHexColor(std::string &out, float r, float g, float b) {
    auto q = [](float c) {
        int v = static_cast<int>(std::lround(std::clamp(c, 0.0f, 1.0f) * 255.0f));
        return std::clamp(v, 0, 255);
    };
    char buf[8];
    std::snprintf(buf, sizeof(buf), "#%02x%02x%02x", q(r), q(g), q(b));
    out.append(buf);
}

struct Vertex {
    float x{0}, y{0};
    float r{0}, g{0}, b{0}, a{1};
};

// One run of primitives that share a style, accumulated as SVG path data so a
// whole wireframe comes out as a handful of elements rather than one per edge.
struct Batch {
    enum class Kind { None, Fill, Stroke } kind{Kind::None};
    float r{0}, g{0}, b{0}, a{1};
    float width{1};
    std::string d;

    bool matches(Kind k, float rr, float gg, float bb, float aa, float w) const {
        if (kind != k) return false;
        const auto same = [](float p, float q) { return std::fabs(p - q) < 1.0f / 512.0f; };
        if (!(same(r, rr) && same(g, gg) && same(b, bb) && same(a, aa))) return false;
        return k != Kind::Stroke || same(width, w);
    }
};

void flush(std::string &out, Batch &batch, std::size_t &elements) {
    if (batch.kind == Batch::Kind::None || batch.d.empty()) {
        batch.kind = Batch::Kind::None;
        batch.d.clear();
        return;
    }
    out += "<path d=\"";
    out += batch.d;
    out += "\" ";
    if (batch.kind == Batch::Kind::Fill) {
        out += "fill=\"";
        appendHexColor(out, batch.r, batch.g, batch.b);
        out += "\" stroke=\"none\"";
        if (batch.a < 0.999f) {
            out += " fill-opacity=\"";
            appendNum(out, batch.a);
            out += "\"";
        }
    } else {
        out += "fill=\"none\" stroke=\"";
        appendHexColor(out, batch.r, batch.g, batch.b);
        out += "\" stroke-width=\"";
        appendNum(out, batch.width);
        out += "\"";
        if (batch.a < 0.999f) {
            out += " stroke-opacity=\"";
            appendNum(out, batch.a);
            out += "\"";
        }
    }
    out += "/>\n";
    ++elements;
    batch.kind = Batch::Kind::None;
    batch.d.clear();
}

} // namespace

// ── theme ────────────────────────────────────────────────────────────────────

void setLightBackground(bool on) { g_light = on; }
bool lightBackground() { return g_light; }

std::array<float, 3> backgroundColor() {
    if (g_light) return {1.0f, 1.0f, 1.0f};
    return {0.1f, 0.1f, 0.12f};
}

void adaptForLight(float &r, float &g, float &b) {
    if (!g_light) return;

    const float mx = std::max({r, g, b});
    const float mn = std::min({r, g, b});

    if (mx - mn < 0.12f) {
        // A neutral: invert its lightness. The range is 0.10 to 0.88 rather
        // than 0 to 1, which keeps the white a wireframe is drawn in off pure
        // black -- softer in print, and no harder than it looked against the
        // grey it was tuned for -- and turns the black the console panel is
        // filled with into a light grey rather than a hole in the page.
        const float v = 1.0f - 0.5f * (mx + mn);
        r = g = b = 0.10f + 0.78f * v;
        return;
    }

    // A hue. Drop half of the darkest channel: against white a colour reads by
    // how far it is *from* white, and a pale hue is nearly nothing.
    r -= 0.5f * mn;
    g -= 0.5f * mn;
    b -= 0.5f * mn;

    // Then cap the luminance, which is what the naturally-bright hues need --
    // a 0.9/0.85/0.1 yellow is invisible on white at full brightness and a
    // perfectly good olive once it is pulled down.
    constexpr float kMaxLum = 0.60f;
    const float lum = 0.2126f * r + 0.7152f * g + 0.0722f * b;
    if (lum > kMaxLum) {
        const float k = kMaxLum / lum;
        r *= k;
        g *= k;
        b *= k;
    }
    r = std::clamp(r, 0.0f, 1.0f);
    g = std::clamp(g, 0.0f, 1.0f);
    b = std::clamp(b, 0.0f, 1.0f);
}

void color3f(float r, float g, float b) {
    adaptForLight(r, g, b);
    glColor3f(r, g, b);
}

void color4f(float r, float g, float b, float a) {
    adaptForLight(r, g, b);
    glColor4f(r, g, b, a);
}

void color3fRaw(float r, float g, float b) { glColor3f(r, g, b); }
void color4fRaw(float r, float g, float b, float a) { glColor4f(r, g, b, a); }

// ── render scale ─────────────────────────────────────────────────────────────

void  setRenderScale(float s) { g_scale = (s > 0.0f) ? s : 1.0f; }
float renderScale() { return g_scale; }
float textScale() { return g_scale; }

void setAxisVisible(bool on) { g_axis = on; }
bool axisVisible() { return g_axis; }

void lineWidth(float w) {
    const float want = w * g_scale;
    if (g_capturing) glPassThrough(kLineWidthMark + want);
    glLineWidth(std::clamp(want, 0.1f, maxGLLineWidth()));
}

void pointSize(float s) {
    const float want = s * g_scale;
    if (g_capturing) glPassThrough(kPointSizeMark + want);
    glPointSize(std::max(0.1f, want));
}

// ── vector capture ───────────────────────────────────────────────────────────

bool capturing() { return g_capturing; }
const CaptureStats &lastCaptureStats() { return g_stats; }

bool beginVectorCapture(std::size_t capacityFloats) {
    if (g_capturing) return false;
    if (capacityFloats < 1024) capacityFloats = 1024;
    g_feedback.assign(capacityFloats, 0.0f);
    glFeedbackBuffer(static_cast<GLsizei>(capacityFloats), GL_3D_COLOR, g_feedback.data());
    glRenderMode(GL_FEEDBACK);
    g_capturing = true;
    g_stats = CaptureStats{};
    return true;
}

CaptureResult buildVectorCapture(std::string &svg, int fbw, int fbh) {
    if (!g_capturing) return CaptureResult::Empty;
    g_capturing = false;

    const GLint used = glRenderMode(GL_RENDER);
    if (used < 0) return CaptureResult::Overflow;
    if (used == 0) return CaptureResult::Empty;
    if (static_cast<std::size_t>(used) > g_feedback.size()) return CaptureResult::Overflow;

    g_stats.floatsUsed = static_cast<std::size_t>(used);

    const int W = std::max(1, fbw);
    const int H = std::max(1, fbh);

    std::string body;
    body.reserve(static_cast<std::size_t>(used) * 4);

    Batch batch;
    std::size_t elements = 0;
    float curLine = 1.0f;
    float curPoint = 1.0f;

    const GLfloat *p   = g_feedback.data();
    const GLfloat *end = p + used;

    // Window coordinates come out of feedback with y up from the bottom-left;
    // SVG has it down from the top-left.
    const auto readVertex = [&](const GLfloat *&q) {
        Vertex v;
        v.x = q[0];
        v.y = static_cast<float>(H) - q[1];
        v.r = q[3];
        v.g = q[4];
        v.b = q[5];
        v.a = q[6];
        q += kVertexFloats;
        return v;
    };

    std::vector<Vertex> poly;

    while (p < end) {
        const int token = static_cast<int>(*p++);

        switch (token) {

        case GL_PASS_THROUGH_TOKEN: {
            if (p >= end) { p = end; break; }
            const float v = *p++;
            if (v >= kPointSizeMark)     curPoint = v - kPointSizeMark;
            else if (v >= kLineWidthMark) curLine  = v - kLineWidthMark;
            break;
        }

        case GL_POINT_TOKEN: {
            if (p + kVertexFloats > end) { p = end; break; }
            const Vertex v = readVertex(p);
            // A point is a disc, not a path segment, so it does not batch.
            flush(body, batch, elements);
            body += "<circle cx=\"";
            appendNum(body, v.x);
            body += "\" cy=\"";
            appendNum(body, v.y);
            body += "\" r=\"";
            appendNum(body, std::max(0.25f, 0.5f * curPoint));
            body += "\" fill=\"";
            appendHexColor(body, v.r, v.g, v.b);
            body += "\"";
            if (v.a < 0.999f) {
                body += " fill-opacity=\"";
                appendNum(body, v.a);
                body += "\"";
            }
            body += "/>\n";
            ++elements;
            ++g_stats.points;
            break;
        }

        case GL_LINE_TOKEN:
        case GL_LINE_RESET_TOKEN: {
            if (p + 2 * kVertexFloats > end) { p = end; break; }
            const Vertex a = readVertex(p);
            const Vertex b = readVertex(p);
            // A gouraud-shaded segment has no single stroke colour; the mean of
            // its ends is the honest answer and, for this viewer, almost always
            // the exact one -- the shading is per-primitive nearly everywhere.
            const float r = 0.5f * (a.r + b.r);
            const float g = 0.5f * (a.g + b.g);
            const float bl = 0.5f * (a.b + b.b);
            const float al = 0.5f * (a.a + b.a);
            if (!batch.matches(Batch::Kind::Stroke, r, g, bl, al, curLine)) {
                flush(body, batch, elements);
                batch.kind = Batch::Kind::Stroke;
                batch.r = r; batch.g = g; batch.b = bl; batch.a = al;
                batch.width = curLine;
            }
            batch.d += 'M';
            appendNum(batch.d, a.x);
            batch.d += ' ';
            appendNum(batch.d, a.y);
            batch.d += 'L';
            appendNum(batch.d, b.x);
            batch.d += ' ';
            appendNum(batch.d, b.y);
            ++g_stats.lines;
            break;
        }

        case GL_POLYGON_TOKEN: {
            if (p >= end) { p = end; break; }
            const int n = static_cast<int>(*p++);
            if (n < 3 || p + static_cast<std::size_t>(n) * kVertexFloats > end) { p = end; break; }
            poly.clear();
            poly.reserve(static_cast<std::size_t>(n));
            float r = 0, g = 0, b = 0, a = 0;
            for (int i = 0; i < n; ++i) {
                const Vertex v = readVertex(p);
                poly.push_back(v);
                r += v.r; g += v.g; b += v.b; a += v.a;
            }
            const float inv = 1.0f / static_cast<float>(n);
            r *= inv; g *= inv; b *= inv; a *= inv;

            if (!batch.matches(Batch::Kind::Fill, r, g, b, a, 0.0f)) {
                flush(body, batch, elements);
                batch.kind = Batch::Kind::Fill;
                batch.r = r; batch.g = g; batch.b = b; batch.a = a;
            }

            // Subpaths of one filled path combine under the nonzero rule, so
            // they have to wind the same way or two adjacent triangles would
            // cancel into a hole. Normalise instead of splitting the batch:
            // the union is what a mesh of abutting triangles should look like,
            // and it comes out without the hairline seams one-path-per-triangle
            // leaves between them.
            double area = 0.0;
            for (int i = 0; i < n; ++i) {
                const Vertex &u = poly[static_cast<std::size_t>(i)];
                const Vertex &w = poly[static_cast<std::size_t>((i + 1) % n)];
                area += static_cast<double>(u.x) * w.y - static_cast<double>(w.x) * u.y;
            }
            if (area < 0.0) std::reverse(poly.begin(), poly.end());

            batch.d += 'M';
            appendNum(batch.d, poly[0].x);
            batch.d += ' ';
            appendNum(batch.d, poly[0].y);
            for (std::size_t i = 1; i < poly.size(); ++i) {
                batch.d += 'L';
                appendNum(batch.d, poly[i].x);
                batch.d += ' ';
                appendNum(batch.d, poly[i].y);
            }
            batch.d += 'Z';
            ++g_stats.polygons;
            break;
        }

        default:
            // An unrecognised token means the stream is no longer being read at
            // a primitive boundary; anything after it would be noise.
            p = end;
            break;
        }
    }
    flush(body, batch, elements);
    g_stats.elements = elements;

    const auto bg = backgroundColor();
    std::ostringstream doc;
    doc << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\"\n"
        << "     width=\"" << W << "\" height=\"" << H << "\"\n"
        << "     viewBox=\"0 0 " << W << " " << H << "\"\n"
        << "     shape-rendering=\"geometricPrecision\">\n"
        << "<g stroke-linecap=\"round\" stroke-linejoin=\"round\">\n";
    {
        std::string rect = "<rect x=\"0\" y=\"0\" width=\"";
        rect += std::to_string(W);
        rect += "\" height=\"";
        rect += std::to_string(H);
        rect += "\" fill=\"";
        appendHexColor(rect, bg[0], bg[1], bg[2]);
        rect += "\"/>\n";
        doc << rect;
    }
    doc << body << "</g>\n</svg>\n";
    svg = doc.str();
    return CaptureResult::Ok;
}

CaptureResult writeVectorCapture(const std::string &path, int fbw, int fbh) {
    std::string svg;
    const CaptureResult r = buildVectorCapture(svg, fbw, fbh);
    if (r != CaptureResult::Ok) return r;

    std::ofstream out(path, std::ios::binary);
    if (!out) return CaptureResult::WriteFailed;
    out << svg;
    out.flush();
    if (!out) return CaptureResult::WriteFailed;
    return CaptureResult::Ok;
}

} // namespace viewer
