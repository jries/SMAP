// Separable image filtering for candidate detection.
//
// The difference-of-Gaussians needs two Gaussians of the same data.  Computing
// them in one traversal halves the memory traffic: the x-pass reads each pixel
// once and writes both convolutions, and the y-pass combines them into the
// final difference without a third pass over the image.
//
// Borders replicate the edge pixel, which is what keeps a bright background
// from producing a spurious step at the image edge.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace smapfit {

// Normalized Gaussian kernel with 2*radius+1 taps.
inline std::vector<float> gauss_kernel(float sigma, int radius) {
    std::vector<float> k(2 * radius + 1);
    double sum = 0.0;
    for (int i = -radius; i <= radius; ++i) {
        const double v = std::exp(-0.5 * double(i) * double(i) / (sigma * sigma));
        k[i + radius] = float(v);
        sum += v;
    }
    for (float& v : k) v = float(v / sum);
    return k;
}

// One-dimensional convolution along x (contiguous) of a single row.
inline void convolve_row(const float* in, float* out, int nx,
                         const float* k, int radius) {
    for (int x = 0; x < nx; ++x) {
        float acc = 0.0f;
        // interior pixels need no clamping; the branch is hoisted out below
        for (int t = -radius; t <= radius; ++t) {
            const int xx = std::min(std::max(x + t, 0), nx - 1);
            acc += in[xx] * k[t + radius];
        }
        out[x] = acc;
    }
}

// Difference of two Gaussians of one frame (ny x nx, row major).
// `tmp` must hold 2 * ny * nx floats.
inline void dog_frame(const float* img, float* out, int ny, int nx,
                      const float* k_narrow, const float* k_wide, int radius,
                      float* tmp) {
    float* ax = tmp;                 // narrow, x-filtered
    float* bx = tmp + ny * nx;       // wide, x-filtered

    // x-pass: both kernels, one traversal of the input
    for (int y = 0; y < ny; ++y) {
        const float* row = img + y * nx;
        float* a = ax + y * nx;
        float* b = bx + y * nx;
        for (int x = 0; x < nx; ++x) {
            float sa = 0.0f, sb = 0.0f;
            const int lo = std::max(x - radius, 0);
            const int hi = std::min(x + radius, nx - 1);
            // clamped tail on the left
            for (int xx = x - radius; xx < lo; ++xx) {
                const int t = xx - x + radius;
                sa += row[0] * k_narrow[t];
                sb += row[0] * k_wide[t];
            }
            for (int xx = lo; xx <= hi; ++xx) {
                const int t = xx - x + radius;
                sa += row[xx] * k_narrow[t];
                sb += row[xx] * k_wide[t];
            }
            // clamped tail on the right
            for (int xx = hi + 1; xx <= x + radius; ++xx) {
                const int t = xx - x + radius;
                sa += row[nx - 1] * k_narrow[t];
                sb += row[nx - 1] * k_wide[t];
            }
            a[x] = sa;
            b[x] = sb;
        }
    }

    // y-pass: accumulate row by row, so every access stays contiguous in x
    for (int y = 0; y < ny; ++y) {
        float* dst = out + y * nx;
        for (int x = 0; x < nx; ++x) dst[x] = 0.0f;
        for (int t = -radius; t <= radius; ++t) {
            const int yy = std::min(std::max(y + t, 0), ny - 1);
            const float wa = k_narrow[t + radius];
            const float wb = k_wide[t + radius];
            const float* a = ax + yy * nx;
            const float* b = bx + yy * nx;
            for (int x = 0; x < nx; ++x) dst[x] += wa * a[x] - wb * b[x];
        }
    }
}

// Difference of two Gaussians as ONE 2-D convolution with the subtracted
// kernel.  The 2-D DoG kernel has rank 2, so it is *not* separable -- but it
// can be applied directly, which trades more arithmetic per pixel ((2r+1)^2
// instead of 4(2r+1)) for a single pass over the image with no temporaries.
// With sigma_wide <= 0 this is a plain 2-D Gaussian (which *is* separable, but
// for a small radius the direct form is still faster).
inline std::vector<float> dog_kernel_2d(float sigma, float sigma_wide,
                                        int radius) {
    const std::vector<float> kn = gauss_kernel(sigma, radius);
    const int n = 2 * radius + 1;
    const std::vector<float> kw =
        sigma_wide > 0.0f ? gauss_kernel(sigma_wide, radius)
                          : std::vector<float>(n, 0.0f);
    std::vector<float> k(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) k[i * n + j] = kn[i] * kn[j] - kw[i] * kw[j];
    return k;
}

inline void dog2d_frame(const float* img, float* out, int ny, int nx,
                        const float* kernel, int radius) {
    const int n = 2 * radius + 1;
    std::vector<float> row(nx);

    for (int y = 0; y < ny; ++y) {
        float* dst = out + y * nx;
        for (int x = 0; x < nx; ++x) dst[x] = 0.0f;

        for (int ty = -radius; ty <= radius; ++ty) {
            const int yy = std::min(std::max(y + ty, 0), ny - 1);
            const float* src = img + yy * nx;
            const float* krow = kernel + (ty + radius) * n;

            for (int tx = -radius; tx <= radius; ++tx) {
                const float w = krow[tx + radius];
                // interior: a plain shifted row, vectorizes cleanly
                const int lo = std::max(-tx, 0);
                const int hi = std::min(nx - tx, nx);
                for (int x = 0; x < lo; ++x) dst[x] += w * src[0];
                for (int x = lo; x < hi; ++x) dst[x] += w * src[x + tx];
                for (int x = hi; x < nx; ++x) dst[x] += w * src[nx - 1];
            }
        }
    }
}

// Single Gaussian of one frame; `tmp` must hold ny * nx floats.
inline void gauss_frame(const float* img, float* out, int ny, int nx,
                        const float* k, int radius, float* tmp) {
    for (int y = 0; y < ny; ++y)
        convolve_row(img + y * nx, tmp + y * nx, nx, k, radius);

    for (int y = 0; y < ny; ++y) {
        float* dst = out + y * nx;
        for (int x = 0; x < nx; ++x) dst[x] = 0.0f;
        for (int t = -radius; t <= radius; ++t) {
            const int yy = std::min(std::max(y + t, 0), ny - 1);
            const float w = k[t + radius];
            const float* src = tmp + yy * nx;
            for (int x = 0; x < nx; ++x) dst[x] += w * src[x];
        }
    }
}

}  // namespace smapfit
