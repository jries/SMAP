// Python bindings for the superresolution renderer.
//
// Separate from _fit3d on purpose: rendering shares no code with fitting, and a
// viewer should not have to carry the fitter.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <stdexcept>

#include "parallel.hpp"
#include "render.hpp"

namespace py = pybind11;

namespace {

using Array = py::array_t<float, py::array::c_style | py::array::forcecast>;

// A per-localization array, or a single value broadcast over all of them.
// Returns the stride to use (1 or 0).
long long broadcast_stride(const Array& a, py::ssize_t n, const char* name) {
    if (a.size() == 1) return 0;
    if (a.size() != n)
        throw std::invalid_argument(std::string(name) +
                                    " must have one value per localization, or one in total");
    return 1;
}

py::tuple render(const Array& x, const Array& y, const Array& sigma_x,
                 const Array& sigma_y, const Array& weight, const Array& color,
                 float x0, float y0, float pixelsize, int nx, int ny,
                 bool gaussian, float roiks, int max_halfwidth, int n_threads) {
    if (x.ndim() != 1 || y.ndim() != 1 || x.size() != y.size())
        throw std::invalid_argument("x and y must be matching 1-D arrays");
    if (nx <= 0 || ny <= 0)
        throw std::invalid_argument("the image must have a positive size");
    if (!(pixelsize > 0))
        throw std::invalid_argument("pixelsize must be positive");

    const py::ssize_t n = x.size();
    const bool colored = color.size() != 0;
    if (colored && (color.ndim() != 2 || color.shape(0) != n || color.shape(1) != 3))
        throw std::invalid_argument("color must have shape (n, 3)");

    Array weight_image({py::ssize_t(ny), py::ssize_t(nx)});
    Array color_image(colored ? std::vector<py::ssize_t>{ny, nx, 3}
                              : std::vector<py::ssize_t>{0, 0, 0});
    std::fill_n(weight_image.mutable_data(), static_cast<size_t>(nx) * ny, 0.0f);
    if (colored)
        std::fill_n(color_image.mutable_data(), static_cast<size_t>(nx) * ny * 3, 0.0f);

    smapfit::RenderInput in{x.data(), y.data(), sigma_x.data(), sigma_y.data(),
                            weight.data(), colored ? color.data() : nullptr, n,
                            gaussian ? broadcast_stride(sigma_x, n, "sigma") : 0,
                            broadcast_stride(weight, n, "weight")};
    if (gaussian && broadcast_stride(sigma_y, n, "sigma") != in.sigma_stride)
        throw std::invalid_argument("sigma_x and sigma_y must broadcast alike");

    smapfit::RenderTarget out{weight_image.mutable_data(),
                              colored ? color_image.mutable_data() : nullptr,
                              nx, ny, x0, y0, 1.0f / pixelsize};

    // One contiguous band of rows per thread; a localization contributes to
    // whichever bands its kernel reaches, and each band writes only its own
    // rows.  With no shared pixels there is nothing to synchronise.
    const int threads = smapfit::resolve_threads(n_threads, ny);
    std::vector<long long> counted(threads, 0);
    {
        py::gil_scoped_release release;
        const int band = (ny + threads - 1) / threads;
        smapfit::parallel_ranges(threads, threads, [&](long long b0, long long b1, int) {
            for (long long b = b0; b < b1; ++b) {
                const int begin = static_cast<int>(b) * band;
                const int end = std::min(begin + band, ny);
                if (begin >= end) continue;
                counted[b] = gaussian
                    ? smapfit::render_gauss_rows(in, out, roiks, max_halfwidth,
                                                 begin, end)
                    : smapfit::render_hist_rows(in, out, begin, end);
            }
        });
    }

    long long total = 0;
    for (long long c : counted) total += c;
    return py::make_tuple(weight_image, colored ? py::object(color_image) : py::none(),
                          total);
}

}  // namespace

PYBIND11_MODULE(_render, m) {
    m.doc() = "Accumulate localizations into a superresolution image.";
    m.def("render", &render, py::arg("x"), py::arg("y"), py::arg("sigma_x"),
          py::arg("sigma_y"), py::arg("weight"), py::arg("color"), py::arg("x0"),
          py::arg("y0"), py::arg("pixelsize"), py::arg("nx"), py::arg("ny"),
          py::arg("gaussian"), py::arg("roiks") = 2.7f,
          py::arg("max_halfwidth") = 512, py::arg("n_threads") = 0,
          "Render localizations; returns (weight, color or None, n_localizations).");
}
