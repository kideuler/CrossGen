#include "SQUADGEN/Sampler.hxx"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

namespace squadgen {
namespace {

// util.spatial_smoothing's lambda and mu (Sec. 1.3, "Taubin weights").
constexpr double kTaubin[2] = {0.4507499669, -0.4720265626};

double orient(const Point &u, const Point &v, const Point &w) {
    return (v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0]);
}

// Whether p-q properly crosses a-b (sq_diffuse.segments_cross).
bool crosses(const Point &p, const Point &q, const Point &a, const Point &b) {
    constexpr double eps = 1e-12;
    return orient(p, q, a) * orient(p, q, b) < -eps && orient(a, b, p) * orient(a, b, q) < -eps;
}

// Dormand-Prince 5(4) as sq_diffuse has it: nodes, stages (the last row is the fifth-order solution, so the seventh
// stage is f at the new point and the next step's first, FSAL) and the fifth- minus fourth-order weights.
constexpr double kC[6] = {1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1.0, 1.0};
constexpr double kA[6][6] = {
    {1.0 / 5},
    {3.0 / 40, 9.0 / 40},
    {44.0 / 45, -56.0 / 15, 32.0 / 9},
    {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729},
    {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656},
    {35.0 / 384, 0.0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84},
};
constexpr double kE[7] = {35.0 / 384 - 5179.0 / 57600, 0.0, 500.0 / 1113 - 7571.0 / 16695, 125.0 / 192 - 393.0 / 640,
                          -2187.0 / 6784 + 92097.0 / 339200, 11.0 / 84 - 187.0 / 2100, -1.0 / 40};

std::vector<float> toFloat(const std::vector<double> &x) { return std::vector<float>(x.begin(), x.end()); }

double rms(const std::vector<double> &v) {
    double s = 0.0;
    for (const double x : v) s += x * x;
    return std::sqrt(s / static_cast<double>(std::max<size_t>(1, v.size())));
}

}  // namespace

int taubinSmooth(std::vector<float> &z, int channels, const std::vector<Point> &tokens,
                 const std::vector<std::array<Point, 2>> &segments, int k, double ringSize, int iterations) {
    const int N = static_cast<int>(tokens.size());
    const size_t C = static_cast<size_t>(channels);
    if (z.size() != static_cast<size_t>(N) * C) throw std::runtime_error("taubinSmooth: the noise is not N x C");
    k = std::min(k, N);
    if (k < 2) return 0;
    // Each token's k nearest, itself first, nearest to farthest (torch.topk, largest = False, sorted).
    std::vector<int> idx(static_cast<size_t>(N) * k);
    std::vector<double> w(static_cast<size_t>(N) * k);
    int blocked = 0;
    std::vector<double> sq(N);
    std::vector<int> order(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const Point d = tokens[j] - tokens[i];
            sq[j] = d[0] * d[0] + d[1] * d[1];
        }
        std::iota(order.begin(), order.end(), 0);
        std::partial_sort(order.begin(), order.begin() + k, order.end(),
                          [&](int a, int b) { return sq[a] < sq[b] || (sq[a] == sq[b] && a < b); });
        const double width = 2.0 * ringSize * ringSize * std::max(sq[order[1]], 1e-20);
        double sum = 0.0;
        for (int m = 0; m < k; ++m) {
            const int j = order[m];
            double wm = std::exp(-sq[j] / width);
            for (const auto &s : segments)
                if (crosses(tokens[i], tokens[j], s[0], s[1])) {
                    wm = 0.0;
                    ++blocked;
                    break;
                }
            idx[static_cast<size_t>(i) * k + m] = j;
            w[static_cast<size_t>(i) * k + m] = wm;
            sum += wm;
        }
        for (int m = 0; m < k; ++m) w[static_cast<size_t>(i) * k + m] /= sum;   // the token itself is never blocked
    }
    std::vector<double> cur(z.begin(), z.end()), next(cur.size());
    for (int it = 0; it < iterations; ++it) {
        const double c = kTaubin[it % 2];
        for (int i = 0; i < N; ++i)
            for (size_t ch = 0; ch < C; ++ch) {
                double avg = 0.0;
                for (int m = 0; m < k; ++m) {
                    const size_t e = static_cast<size_t>(i) * k + m;
                    avg += w[e] * cur[static_cast<size_t>(idx[e]) * C + ch];
                }
                const size_t o = static_cast<size_t>(i) * C + ch;
                next[o] = c * avg + (1.0 - c) * cur[o];
            }
        cur.swap(next);
    }
    for (size_t o = 0; o < z.size(); ++o) z[o] = static_cast<float>(cur[o]);
    return blocked;
}

Sampler::Result Sampler::integrate(const Velocity &f, const std::vector<float> &x0, const Options &opts,
                                   const Observer &observe) {
    const size_t M = x0.size();
    Result res;
    std::vector<double> x(x0.begin(), x0.end());
    std::vector<float> xf = x0, vf;
    auto eval = [&](const std::vector<double> &y, double t, std::vector<double> &k) {
        xf.assign(y.begin(), y.end());
        f(xf, t, vf);
        if (vf.size() != M) throw std::runtime_error("SQuadGen: the velocity has the wrong size");
        k.assign(vf.begin(), vf.end());
        ++res.evaluations;
    };
    auto show = [&](double t, double h, const std::vector<double> &v) {
        if (!observe) return true;
        const std::vector<float> xs = toFloat(x), vs = toFloat(v);
        Progress p;
        p.t = t;
        p.h = h;
        p.steps = res.steps;
        p.rejected = res.rejected;
        p.evaluations = res.evaluations;
        p.x = &xs;
        p.v = &vs;
        return observe(p);
    };
    auto finish = [&](double t, bool done) {
        res.x = toFloat(x);
        res.t = t;
        res.finished = done;
        return res;
    };

    if (opts.method != Method::Dopri5) {
        const int steps = std::max(1, opts.steps);
        const double h = 1.0 / steps;
        std::vector<double> v, v2, y(M);
        for (int i = 0; i < steps; ++i) {
            eval(x, i * h, v);
            if (i == 0 && !show(0.0, h, v)) return finish(0.0, false);
            if (opts.method == Method::Heun) {
                for (size_t e = 0; e < M; ++e) y[e] = x[e] + h * v[e];
                eval(y, (i + 1) * h, v2);
                for (size_t e = 0; e < M; ++e) v[e] = 0.5 * (v[e] + v2[e]);
            }
            for (size_t e = 0; e < M; ++e) x[e] += h * v[e];
            ++res.steps;
            const double t = (i + 1 == steps) ? 1.0 : (i + 1) * h;
            if (!show(t, h, v)) return finish(t, i + 1 == steps);
        }
        return finish(1.0, true);
    }

    // Dopri5 (sq_diffuse.dopri5).
    const double t0 = 0.0, t1 = 1.0, atol = opts.atol, rtol = opts.rtol;
    std::vector<double> k1, kk, y(M), scaled(M);
    eval(x, t0, k1);
    for (size_t e = 0; e < M; ++e) scaled[e] = x[e] / (atol + std::fabs(x[e]) * rtol);
    const double d0 = rms(scaled);
    for (size_t e = 0; e < M; ++e) scaled[e] = k1[e] / (atol + std::fabs(x[e]) * rtol);
    const double d1 = rms(scaled);
    const double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
    for (size_t e = 0; e < M; ++e) y[e] = x[e] + h0 * k1[e];
    eval(y, t0 + h0, kk);
    for (size_t e = 0; e < M; ++e) scaled[e] = (kk[e] - k1[e]) / (atol + std::fabs(x[e]) * rtol);
    const double d2 = rms(scaled) / h0;
    const double dm = std::max(d1, d2);
    double h = std::min(100.0 * h0, dm <= 1e-15 ? std::max(1e-6, h0 * 1e-3) : std::pow(0.01 / dm, 0.2));
    double t = t0;
    if (!show(t, h, k1)) return finish(t, false);

    std::vector<std::vector<double>> ks(7, std::vector<double>(M));
    std::vector<double> err(M);
    for (int it = 0; it < opts.maxSteps; ++it) {
        if (t1 - t <= 1e-12) return finish(1.0, true);
        h = std::min(h, t1 - t);
        ks[0] = k1;
        for (int s = 0; s < 6; ++s) {
            for (size_t e = 0; e < M; ++e) {
                double acc = 0.0;
                for (int j = 0; j <= s; ++j)
                    if (kA[s][j] != 0.0) acc += kA[s][j] * ks[j][e];
                y[e] = x[e] + h * acc;
            }
            eval(y, t + kC[s] * h, ks[s + 1]);
        }
        for (size_t e = 0; e < M; ++e) {
            double acc = 0.0;
            for (int j = 0; j < 7; ++j)
                if (kE[j] != 0.0) acc += kE[j] * ks[j][e];
            err[e] = h * acc / (atol + rtol * std::max(std::fabs(x[e]), std::fabs(y[e])));
        }
        const double ratio = rms(err);
        const double used = h;
        if (ratio <= 1.0) {
            t += h;
            x.swap(y);
            k1 = ks[6];
            ++res.steps;
        } else {
            ++res.rejected;
        }
        h *= ratio == 0.0 ? 10.0 : std::min(10.0, std::max(0.9 * std::pow(ratio, -0.2), ratio > 1.0 ? 0.2 : 1.0));
        if (ratio <= 1.0 && !show(t, used, k1)) return finish(t, t1 - t <= 1e-12);
    }
    throw std::runtime_error("dopri5: no convergence in " + std::to_string(opts.maxSteps) + " steps (t = " +
                             std::to_string(t) + ", h = " + std::to_string(h) + ")");
}

const char *Sampler::name(Method m) {
    switch (m) {
        case Method::Dopri5: return "dopri5";
        case Method::Heun:   return "heun";
        case Method::Euler:  return "euler";
    }
    return "?";
}

}  // namespace squadgen
