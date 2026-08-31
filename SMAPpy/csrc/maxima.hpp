// Strict local-maximum search, following SMAP's maximumfind2.c.
//
// A pixel is a maximum when it is strictly greater than all eight neighbours,
// so plateaus yield nothing.  The single pass with short-circuit evaluation
// touches only one or two neighbours for the vast majority of pixels, which is
// what makes it faster than a filter-based search: those have to compute a
// neighbourhood maximum for every pixel before comparing.
#pragma once

#include <cstdint>
#include <vector>

namespace smappy {

struct Maximum {
    int32_t frame, y, x;
    float value;
};

// Append all strict maxima of one frame (row-major, ny x nx) to `out`.
// The one-pixel border is skipped, as in the original.
inline void find_maxima_frame(const float* img, int ny, int nx, int frame,
                              float threshold, std::vector<Maximum>& out) {
    for (int y = 1; y < ny - 1; ++y) {
        const float* row = img + y * nx;
        const float* above = row - nx;
        const float* below = row + nx;
        for (int x = 1; x < nx - 1; ++x) {
            const float v = row[x];
            if (v <= threshold) continue;
            if (v > row[x - 1] && v > row[x + 1] &&
                v > above[x - 1] && v > above[x] && v > above[x + 1] &&
                v > below[x - 1] && v > below[x] && v > below[x + 1]) {
                out.push_back({frame, y, x, v});
                ++x;  // the next pixel is smaller, so it cannot be a maximum
            }
        }
    }
}

}  // namespace smappy
