// PSF models for the maximum-likelihood fitter.
//
// Each model provides, for one ROI:
//   NV                       number of fit parameters
//   init(data, sz, theta, maxjump)      starting values
//   prepare(theta, sz)                  per-iteration precomputation
//   value(ix, iy, theta, dudt, &model)  model value and its derivatives
//   clamp(theta, sz)                    keep parameters in a sane range
//
// Parameter order is (x, y, photons, background, [z | sigma | sigma_x, sigma_y]).
// `ix` is the ROI column and `iy` the row: the data of one ROI is row-major,
// element (iy, ix) at data[iy*sz + ix], and theta[0] is the x coordinate.
//
// The model formulae are transcribed from the SMAP / GPUmleFit_LM sources
// (CPUgaussLib.cpp, CPUsplineLib.cpp, CPUmleFit_LM.cpp).  The one deliberate
// change is that x and y are no longer restricted to the central half of the
// ROI: the multi-channel fitter dropped that restriction and we follow it, so
// poorly centred candidates produce a visibly bad fit instead of one silently
// parked on the clamp boundary.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>

namespace smappy {

constexpr float PI_F = 3.141592f;

// ---------------------------------------------------------------- Gaussian PSF
inline float int_gauss_1d(int ii, float x, float sigma) {
    const float norm = 1.0f / 2.0f / sigma / sigma;
    return 0.5f * (std::erf((ii - x + 0.5f) * std::sqrt(norm)) -
                   std::erf((ii - x - 0.5f) * std::sqrt(norm)));
}

inline void d_int_gauss_1d(int ii, float x, float sigma, float N, float PSFy,
                           float* dudt) {
    const float a = std::exp(-0.5f * (ii + 0.5f - x) * (ii + 0.5f - x) / (sigma * sigma));
    const float b = std::exp(-0.5f * (ii - 0.5f - x) * (ii - 0.5f - x) / (sigma * sigma));
    *dudt = -N / std::sqrt(2.0f * PI_F) / sigma * (a - b) * PSFy;
}

inline void d_int_gauss_1d_sigma(int ii, float x, float sx, float N, float PSFy,
                                 float* dudt) {
    const float ax = std::exp(-0.5f * (ii + 0.5f - x) * (ii + 0.5f - x) / (sx * sx));
    const float bx = std::exp(-0.5f * (ii - 0.5f - x) * (ii - 0.5f - x) / (sx * sx));
    *dudt = -N / std::sqrt(2.0f * PI_F) / sx / sx *
            (ax * (ii - x + 0.5f) - bx * (ii - x - 0.5f)) * PSFy;
}

// Centre of mass and min/max of a ROI, used for the starting values.
inline void center_of_mass(const float* data, int sz, float* x, float* y) {
    float sx = 0.0f, sy = 0.0f, sum = 0.0f;
    for (int iy = 0; iy < sz; ++iy)
        for (int ix = 0; ix < sz; ++ix) {
            const float v = data[iy * sz + ix];
            sx += v * ix;
            sy += v * iy;
            sum += v;
        }
    *x = sx / sum;
    *y = sy / sum;
}

inline void max_min(const float* data, int sz, float* maxn, float* minbg) {
    *maxn = 0.0f;
    *minbg = 1e10f;
    for (int i = 0; i < sz * sz; ++i) {
        *maxn = std::max(*maxn, data[i]);
        *minbg = std::min(*minbg, data[i]);
    }
}

// Gaussian with one free width.  theta = (x, y, N, bg, sigma)
struct GaussFree {
    static constexpr int NV = 5;
    float sigma_start;

    explicit GaussFree(float sigma = 1.0f) : sigma_start(sigma) {}

    void init(const float* data, int sz, float* theta, float* maxjump) const {
        float maxn, minbg;
        center_of_mass(data, sz, &theta[0], &theta[1]);
        max_min(data, sz, &maxn, &minbg);
        theta[3] = std::max(minbg, 0.01f);
        theta[2] = std::max(0.0f, (maxn - theta[3]) * 2 * PI_F * sigma_start * sigma_start);
        theta[4] = sigma_start;

        maxjump[0] = 1.0f; maxjump[1] = 1.0f;
        maxjump[2] = std::max(theta[2], 100.0f);
        maxjump[3] = std::max(theta[3], 20.0f);
        maxjump[4] = 0.5f;
    }

    void prepare(const float*, int) const {}

    void value(int ix, int iy, const float* theta, float* dudt, float* model) const {
        const float px = int_gauss_1d(ix, theta[0], theta[4]);
        const float py = int_gauss_1d(iy, theta[1], theta[4]);
        *model = theta[3] + theta[2] * px * py;
        d_int_gauss_1d(ix, theta[0], theta[4], theta[2], py, &dudt[0]);
        d_int_gauss_1d(iy, theta[1], theta[4], theta[2], px, &dudt[1]);
        float dsx, dsy;
        d_int_gauss_1d_sigma(ix, theta[0], theta[4], theta[2], py, &dsx);
        d_int_gauss_1d_sigma(iy, theta[1], theta[4], theta[2], px, &dsy);
        dudt[4] = dsx + dsy;
        dudt[2] = px * py;
        dudt[3] = 1.0f;
    }

    void clamp(float* theta, int sz) const {
        theta[2] = std::max(theta[2], 1.0f);
        theta[3] = std::max(theta[3], 0.01f);
        theta[4] = std::max(theta[4], 0.0f);
        theta[4] = std::min(theta[4], sz / 2.0f);
    }
};

// Elliptical Gaussian.  theta = (x, y, N, bg, sigma_x, sigma_y)
struct GaussXY {
    static constexpr int NV = 6;
    float sigma_start;

    explicit GaussXY(float sigma = 1.0f) : sigma_start(sigma) {}

    void init(const float* data, int sz, float* theta, float* maxjump) const {
        float maxn, minbg;
        center_of_mass(data, sz, &theta[0], &theta[1]);
        max_min(data, sz, &maxn, &minbg);
        theta[3] = std::max(minbg, 0.01f);
        theta[2] = std::max(0.0f, (maxn - theta[3]) * 2 * PI_F * sigma_start * sigma_start);
        theta[4] = sigma_start;
        theta[5] = sigma_start;

        maxjump[0] = 1.0f; maxjump[1] = 1.0f;
        maxjump[2] = std::max(theta[2], 100.0f);
        maxjump[3] = std::max(theta[3], 20.0f);
        maxjump[4] = 0.5f; maxjump[5] = 0.5f;
    }

    void prepare(const float*, int) const {}

    void value(int ix, int iy, const float* theta, float* dudt, float* model) const {
        const float px = int_gauss_1d(ix, theta[0], theta[4]);
        const float py = int_gauss_1d(iy, theta[1], theta[5]);
        *model = theta[3] + theta[2] * px * py;
        d_int_gauss_1d(ix, theta[0], theta[4], theta[2], py, &dudt[0]);
        d_int_gauss_1d(iy, theta[1], theta[5], theta[2], px, &dudt[1]);
        d_int_gauss_1d_sigma(ix, theta[0], theta[4], theta[2], py, &dudt[4]);
        d_int_gauss_1d_sigma(iy, theta[1], theta[5], theta[2], px, &dudt[5]);
        dudt[2] = px * py;
        dudt[3] = 1.0f;
    }

    void clamp(float* theta, int sz) const {
        theta[2] = std::max(theta[2], 1.0f);
        theta[3] = std::max(theta[3], 0.01f);
        theta[4] = std::max(theta[4], sigma_start / 10.0f);
        theta[5] = std::max(theta[5], sigma_start / 10.0f);
    }
};

// ------------------------------------------------------------------ cspline PSF
// The 64 monomials of the tricubic polynomial, ordered i = 16*pz + 4*py + px.
inline void compute_delta3d(float dx, float dy, float dz, float* f, float* dfx,
                            float* dfy, float* dfz) {
    std::memset(f, 0, 64 * sizeof(float));
    std::memset(dfx, 0, 64 * sizeof(float));
    std::memset(dfy, 0, 64 * sizeof(float));
    std::memset(dfz, 0, 64 * sizeof(float));

    float cz = 1.0f;
    for (int i = 0; i < 4; ++i) {
        float cy = 1.0f;
        for (int j = 0; j < 4; ++j) {
            float cx = 1.0f;
            for (int k = 0; k < 4; ++k) {
                f[i * 16 + j * 4 + k] = cz * cy * cx;
                if (k < 3) dfx[i * 16 + j * 4 + k + 1] = (k + 1) * cz * cy * cx;
                if (j < 3) dfy[i * 16 + (j + 1) * 4 + k] = (j + 1) * cz * cy * cx;
                if (i < 3) dfz[(i + 1) * 16 + j * 4 + k] = (i + 1) * cz * cy * cx;
                cx *= dx;
            }
            cy *= dy;
        }
        cz *= dz;
    }
}

// Cubic spline PSF.  theta = (x, y, N, bg, z)
// Coefficients are (64, nz, ny, nx) C-contiguous, x fastest; see
// smappy.io.calibration for how the SMAP calibration is converted.
struct CSpline {
    static constexpr int NV = 5;

    const float* coeff;
    int nx, ny, nz;
    float z_start;

    // per-fit state, refreshed by prepare()
    mutable float delta_f[64], delta_dx[64], delta_dy[64], delta_dz[64];
    mutable int xstart, ystart, zstart, off_x, off_y;

    CSpline(const float* c, int nx_, int ny_, int nz_, float z0)
        : coeff(c), nx(nx_), ny(ny_), nz(nz_), z_start(z0) {}

    void init(const float* data, int sz, float* theta, float* maxjump) const {
        float maxn, minbg;
        center_of_mass(data, sz, &theta[0], &theta[1]);
        max_min(data, sz, &maxn, &minbg);
        theta[3] = std::max(minbg, 0.01f);
        const float centre = coeff[(((nz / 2) * ny) + ny / 2) * nx + nx / 2];
        theta[2] = (maxn - theta[3]) / centre * 4.0f;
        theta[4] = z_start;

        maxjump[0] = 1.0f; maxjump[1] = 1.0f;
        maxjump[2] = std::max(theta[2], 100.0f);
        maxjump[3] = std::max(theta[3], 20.0f);
        maxjump[4] = std::max(nz / 3.0f, 2.0f);
    }

    void prepare(const float* theta, int sz) const {
        float xc = -1.0f * ((theta[0] - sz / 2.0f) + 0.5f);
        float yc = -1.0f * ((theta[1] - sz / 2.0f) + 0.5f);
        off_x = static_cast<int>(std::floor((nx + 1.0f - sz) / 2.0f));
        off_y = static_cast<int>(std::floor((ny + 1.0f - sz) / 2.0f));
        xstart = static_cast<int>(std::floor(xc));
        ystart = static_cast<int>(std::floor(yc));
        zstart = static_cast<int>(std::floor(theta[4]));
        compute_delta3d(xc - xstart, yc - ystart, theta[4] - zstart, delta_f,
                        delta_dx, delta_dy, delta_dz);
    }

    void value(int ix, int iy, const float* theta, float* dudt, float* model) const {
        int xc = std::min(std::max(ix + xstart + off_x, 0), nx - 1);
        int yc = std::min(std::max(iy + ystart + off_y, 0), ny - 1);
        int zc = std::min(std::max(zstart, 0), nz - 1);

        const int base = (zc * ny + yc) * nx + xc;
        const int stride = nx * ny * nz;
        float temp = 0.0f, dx = 0.0f, dy = 0.0f, dz = 0.0f;
        for (int i = 0; i < 64; ++i) {
            const float c = coeff[i * stride + base];
            temp += delta_f[i] * c;
            dx += delta_dx[i] * c;
            dy += delta_dy[i] * c;
            dz += delta_dz[i] * c;
        }
        dudt[0] = -theta[2] * dx;
        dudt[1] = -theta[2] * dy;
        dudt[4] = theta[2] * dz;
        dudt[2] = temp;
        dudt[3] = 1.0f;
        *model = theta[3] + theta[2] * temp;
    }

    void clamp(float* theta, int) const {
        theta[2] = std::max(theta[2], 1.0f);
        theta[3] = std::max(theta[3], 0.01f);
        theta[4] = std::max(theta[4], 0.0f);
        theta[4] = std::min(theta[4], static_cast<float>(nz));
    }
};

}  // namespace smappy
