#ifndef __SQUADGEN_SAMPLER_HXX__
#define __SQUADGEN_SAMPLER_HXX__

#include <array>
#include <functional>
#include <vector>

#include "mesh/Mesh.hxx"

namespace squadgen {

// Appendix C (sq_diffuse.taubin_smooth): Taubin smoothing of the initial noise z, N x C row-major, over each token's k
// nearest tokens (itself included), with Gaussian weights whose width is ringSize times the distance to its nearest
// other token, and ten alternating lambda / mu passes (util.spatial_smoothing's 0.4507499669 / -0.4720265626). In
// place of the release's normal term, a pair whose segment properly crosses dS gets no weight, so that noise is not
// correlated across a slot or a thin wall. Returns how many of the N k pairs were blocked so.
int taubinSmooth(std::vector<float> &z, int channels, const std::vector<Point> &tokens,
                 const std::vector<std::array<Point, 2>> &segments, int k = 32, double ringSize = 2.0,
                 int iterations = 10);

// The deterministic ODE SQ-Diffuse samples by (Sec. 3.4, the release's Sampler.sample_ode): dx/dt = v(x, t) from the
// noise at t = 0 to the latent at t = 1, on normalised latents.
//
//   Dopri5  Dormand-Prince 5(4) with torchdiffeq's initial step, error norm (the RMS over every element of
//           err / (atol + rtol max(|x|, |x_new|))) and controller (safety 0.9, growth at most 10, shrink at least
//           0.2, exponent 1/5): sq_diffuse.dopri5, line for line, the last step cut to land on t = 1. The release's
//           inference, at atol 1e-6 and rtol 1e-3: 100-150 evaluations a sample.
//   Heun    the trapezoidal predictor-corrector on `steps` equal steps, two evaluations each (the trainer's
//           evaluation sampler)
//   Euler   one evaluation a step
//
// The state is kept in double; each evaluation is handed the float32 the networks take.
class Sampler {
public:
    enum class Method { Dopri5, Heun, Euler };

    struct Options {
        Method method = Method::Dopri5;
        int steps = 12;            // Heun / Euler
        double atol = 1e-6;        // Dopri5
        double rtol = 1e-3;
        int maxSteps = 10000;
    };

    // What an observer is shown: once at t = 0, with the first evaluation, and then after every accepted step. `v` is
    // the velocity at (x, t) for Dopri5 (its FSAL stage), the step's own average for Heun and the start's for Euler,
    // so x + (1 - t) v is the latent the flow predicts at t = 1 -- what a live picture of the solve decodes.
    struct Progress {
        double t = 0.0, h = 0.0;
        int steps = 0, rejected = 0, evaluations = 0;
        const std::vector<float> *x = nullptr;
        const std::vector<float> *v = nullptr;
    };
    // Returning false stops the solve where it is.
    using Observer = std::function<bool(const Progress &)>;
    // One evaluation: v = f(x, t), resized by the callee.
    using Velocity = std::function<void(const std::vector<float> &x, double t, std::vector<float> &v)>;

    struct Result {
        std::vector<float> x;      // at t = 1, or where an observer stopped it
        double t = 0.0;
        int evaluations = 0, steps = 0, rejected = 0;
        bool finished = false;
    };

    static Result integrate(const Velocity &f, const std::vector<float> &x0, const Options &opts,
                            const Observer &observe = {});

    static const char *name(Method m);
};

}  // namespace squadgen

#endif // __SQUADGEN_SAMPLER_HXX__
