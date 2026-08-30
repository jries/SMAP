// Python binding for the grouper's linking core.  See group.hpp.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <stdexcept>

#include "group.hpp"

namespace py = pybind11;

namespace {

using Doubles = py::array_t<double, py::array::c_style | py::array::forcecast>;
using Frames = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;

py::tuple connect(const Doubles& x, const Doubles& y, const Frames& frame,
                  double dx, int64_t dt) {
    if (x.ndim() != 1 || y.size() != x.size() || frame.size() != x.size())
        throw std::invalid_argument("x, y and frame must be matching 1-D arrays");

    const py::ssize_t n = x.size();
    py::array_t<int64_t> list(n);
    std::fill_n(list.mutable_data(), n, int64_t(0));

    int64_t groups = 0;
    {
        py::gil_scoped_release release;   // sequential, but long-running
        groups = smapfit::connect_single(x.data(), y.data(), frame.data(), n, dx,
                                         dt, list.mutable_data());
    }
    return py::make_tuple(list, groups);
}

}  // namespace

PYBIND11_MODULE(_group, m) {
    m.doc() = "Frame-to-frame linking of localizations.";
    m.def("connect", &connect, py::arg("x"), py::arg("y"), py::arg("frame"),
          py::arg("dx"), py::arg("dt"),
          "Link one block sorted by (frame, x); returns (group_ids, n_groups).");
}
