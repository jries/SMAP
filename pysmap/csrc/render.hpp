// Accumulating localizations into a superresolution image.
//
// One kernel routine covers every Gaussian mode: the caller supplies a sigma
// per localization (or one for all), and the histogram is the sigma -> 0 limit.
//
// Two departures from SMAP's gaussrender, both because they make the result
// exact rather than merely fast:
//
//  * the kernel is the **pixel integral** of the Gaussian, erf(right) -
//    erf(left), not a point sample of it.  A rendered pixel is an integral over
//    its area, and sampling breaks down exactly where sigma approaches one
//    pixel -- which is where the renderer spends most of its time.  It is also
//    what lets the histogram be the limiting case instead of a separate model.
//  * the kernel is normalised by its own sum, not by an erf correction for the
//    truncation.  Every localization then contributes exactly N to the image,
//    whatever its sigma or its subpixel position, so the Gaussian and histogram
//    modes carry the same total intensity and are directly comparable.
//
// The kernel is separable, so a ROI of (2d+1)^2 pixels costs 2(2d+2) erf calls
// and (2d+1)^2 multiply-adds.  SMAP instead looks up a 601x601 template with
// nearest-neighbour sampling; that trades the quantisation error for nothing we
// need, since the outer product dominates the cost either way.
//
// Threading: a thread owns a contiguous range of output rows and simply skips
// what falls outside it, so no two threads ever write the same pixel.  No
// locks, no per-thread image copies.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace smapfit {

// Localizations, already in the coordinate system the caller wants; the render
// converts to pixels itself so no temporary arrays are needed.
struct RenderInput {
    const float* x;
    const float* y;
    const float* sigma_x;   // may be a single value, see stride
    const float* sigma_y;
    const float* weight;    // intensity per localization; may be a single value
    const float* color;     // (n, 3) or null
    long long n;
    long long sigma_stride; // 1 for per-localization, 0 for one for all
    long long weight_stride;
};

struct RenderTarget {
    float* weight;          // (ny, nx)
    float* color;           // (ny, nx, 3) or null
    int nx, ny;
    float x0, y0, inv_pixelsize;
};

// The pixel-integrated Gaussian over [-d, d] around the centre pixel, written
// to k[0..2d]; returns its sum.  `offset` is the position of the localization
// relative to the centre of the centre pixel, in pixels.
inline float gauss_kernel_1d(float* k, int d, float offset, float sigma) {
    const float scale = 1.0f / (1.41421356237f * sigma);
    float prev = std::erf((-d - 0.5f - offset) * scale);
    float sum = 0.0f;
    for (int i = 0; i <= 2 * d; ++i) {
        const float next = std::erf((i - d + 0.5f - offset) * scale);
        k[i] = 0.5f * (next - prev);
        sum += k[i];
        prev = next;
    }
    return sum;
}

// Half-width of the ROI in pixels, following SMAP's roiks * sigma + 1.
inline int kernel_halfwidth(float sigma, float roiks, int limit) {
    const int d = static_cast<int>(roiks * sigma + 1.0f);
    return std::min(std::max(d, 0), limit);
}

// Accumulate into rows [row_begin, row_end) of the target.  Returns the number
// of localizations whose centre pixel lies in those rows (and in the image),
// so that summing over threads counts every visible localization once.
inline long long render_gauss_rows(const RenderInput& in, const RenderTarget& out,
                                   float roiks, int max_halfwidth,
                                   int row_begin, int row_end) {
    std::vector<float> kx, ky;
    long long counted = 0;

    for (long long i = 0; i < in.n; ++i) {
        const float px = (in.x[i] - out.x0) * out.inv_pixelsize;
        const float py = (in.y[i] - out.y0) * out.inv_pixelsize;
        if (!(std::isfinite(px) && std::isfinite(py))) continue;

        const int cx = static_cast<int>(std::floor(px));
        const int cy = static_cast<int>(std::floor(py));
        const float sx = in.sigma_x[i * in.sigma_stride] * out.inv_pixelsize;
        const float sy = in.sigma_y[i * in.sigma_stride] * out.inv_pixelsize;

        // sigma below a hundredth of a pixel is a histogram entry: the pixel
        // integral has collapsed onto one pixel anyway
        const int dx = sx > 0.01f ? kernel_halfwidth(sx, roiks, max_halfwidth) : 0;
        const int dy = sy > 0.01f ? kernel_halfwidth(sy, roiks, max_halfwidth) : 0;

        const int col_lo = std::max(cx - dx, 0), col_hi = std::min(cx + dx, out.nx - 1);
        const int row_lo = std::max(std::max(cy - dy, 0), row_begin);
        const int row_hi = std::min(std::min(cy + dy, out.ny - 1), row_end - 1);
        if (cy >= row_begin && cy < row_end && cy >= 0 && cy < out.ny &&
            cx >= 0 && cx < out.nx)
            ++counted;
        if (col_lo > col_hi || row_lo > row_hi) continue;

        // The kernel is built over the full ROI even where the image clips it,
        // so a localization at the border loses the part that falls outside
        // instead of having it redistributed inwards -- and so that the result
        // does not depend on where the thread boundaries are.
        float norm = in.weight[i * in.weight_stride];
        if (dx > 0) {
            kx.resize(2 * dx + 1);
            norm /= gauss_kernel_1d(kx.data(), dx, px - cx - 0.5f, sx);
        }
        if (dy > 0) {
            ky.resize(2 * dy + 1);
            norm /= gauss_kernel_1d(ky.data(), dy, py - cy - 0.5f, sy);
        }

        const float* cols = in.color ? in.color + 3 * i : nullptr;
        for (int row = row_lo; row <= row_hi; ++row) {
            const float wy = dy > 0 ? ky[row - cy + dy] * norm : norm;
            float* wrow = out.weight + static_cast<long long>(row) * out.nx;
            float* crow = out.color
                ? out.color + 3 * static_cast<long long>(row) * out.nx : nullptr;
            for (int col = col_lo; col <= col_hi; ++col) {
                const float w = dx > 0 ? wy * kx[col - cx + dx] : wy;
                wrow[col] += w;
                if (crow) {
                    float* p = crow + 3 * col;
                    p[0] += w * cols[0];
                    p[1] += w * cols[1];
                    p[2] += w * cols[2];
                }
            }
        }
    }
    return counted;
}

// The histogram: every localization lands in one pixel.  Kept separate because
// it needs no kernel at all, and it is the mode used on the densest data.
inline long long render_hist_rows(const RenderInput& in, const RenderTarget& out,
                                  int row_begin, int row_end) {
    long long counted = 0;
    for (long long i = 0; i < in.n; ++i) {
        const float px = (in.x[i] - out.x0) * out.inv_pixelsize;
        const float py = (in.y[i] - out.y0) * out.inv_pixelsize;
        if (!(std::isfinite(px) && std::isfinite(py))) continue;
        const int col = static_cast<int>(std::floor(px));
        const int row = static_cast<int>(std::floor(py));
        if (col < 0 || col >= out.nx || row < row_begin || row >= row_end ||
            row < 0 || row >= out.ny)
            continue;
        ++counted;
        const long long p = static_cast<long long>(row) * out.nx + col;
        const float w = in.weight[i * in.weight_stride];
        out.weight[p] += w;
        if (out.color) {
            const float* c = in.color + 3 * i;
            out.color[3 * p + 0] += w * c[0];
            out.color[3 * p + 1] += w * c[1];
            out.color[3 * p + 2] += w * c[2];
        }
    }
    return counted;
}

}  // namespace smapfit
