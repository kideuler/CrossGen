#include "python/PyUtil.hxx"

#include <cstdint>
#include <cstring>

namespace pycg {

std::mutex &computeMutex() {
    static std::mutex m;
    return m;
}

namespace {

// numpy, imported once. Borrowed from this cache; null with an ImportError set
// if the interpreter has none.
PyObject *numpy() {
    static PyObject *np = nullptr;
    if (!np) {
        np = PyImport_ImportModule("numpy");
        if (!np) {
            PyErr_SetString(PyExc_ImportError,
                            "crossgen returns numpy arrays, and numpy could not be imported");
            return nullptr;
        }
    }
    return np;
}

// numpy.frombuffer over a bytearray holding a copy of `bytes` bytes of `data`,
// reshaped to (rows, cols), or (rows,) when cols is 0. The bytearray is what
// makes the array writeable and its owner.
PyObject *makeArray(const void *data, size_t bytes, const char *dtype, Py_ssize_t rows,
                    Py_ssize_t cols) {
    PyObject *np = numpy();
    if (!np) return nullptr;
    PyObject *buf = PyByteArray_FromStringAndSize(nullptr, static_cast<Py_ssize_t>(bytes));
    if (!buf) return nullptr;
    if (bytes > 0) std::memcpy(PyByteArray_AS_STRING(buf), data, bytes);
    PyObject *flat = PyObject_CallMethod(np, "frombuffer", "Os", buf, dtype);
    Py_DECREF(buf);
    if (!flat || cols <= 0) return flat;
    PyObject *shaped = PyObject_CallMethod(flat, "reshape", "nn", rows, cols);
    Py_DECREF(flat);
    return shaped;
}

// numpy.ascontiguousarray(o, dtype) with ndim and the column count checked,
// and its buffer. The caller releases `view` and `arr`.
PyObject *contiguous(PyObject *o, const char *dtype, int cols, const char *what, Py_buffer &view) {
    PyObject *np = numpy();
    if (!np) return nullptr;
    PyObject *arr = PyObject_CallMethod(np, "ascontiguousarray", "Os", o, dtype);
    if (!arr) return nullptr;
    if (PyObject_GetBuffer(arr, &view, PyBUF_C_CONTIGUOUS) < 0) {
        Py_DECREF(arr);
        return nullptr;
    }
    const int ndim = cols > 0 ? 2 : 1;
    const bool shapeOk = view.ndim == ndim && (cols <= 0 || view.shape[1] == cols) &&
                         view.itemsize == 8;
    if (!shapeOk) {
        if (cols > 0)
            PyErr_Format(PyExc_ValueError, "%s must have shape (n, %d)", what, cols);
        else
            PyErr_Format(PyExc_ValueError, "%s must be one-dimensional", what);
        PyBuffer_Release(&view);
        Py_DECREF(arr);
        return nullptr;
    }
    return arr;
}

}  // namespace

PyObject *doubleArray(const double *data, Py_ssize_t rows, Py_ssize_t cols) {
    const size_t n = static_cast<size_t>(rows) * static_cast<size_t>(cols > 0 ? cols : 1);
    return makeArray(data, n * sizeof(double), "float64", rows, cols);
}

PyObject *indexArray(const int *data, Py_ssize_t rows, Py_ssize_t cols) {
    const size_t n = static_cast<size_t>(rows) * static_cast<size_t>(cols > 0 ? cols : 1);
    std::vector<std::int64_t> wide(n);
    for (size_t i = 0; i < n; ++i) wide[i] = data[i];
    return makeArray(wide.data(), n * sizeof(std::int64_t), "int64", rows, cols);
}

PyObject *doubleArray(const std::vector<double> &v) {
    return doubleArray(v.data(), static_cast<Py_ssize_t>(v.size()), 0);
}

PyObject *indexArray(const std::vector<int> &v) {
    return indexArray(v.data(), static_cast<Py_ssize_t>(v.size()), 0);
}

PyObject *pointArray(const std::vector<Point> &points) {
    return doubleArray(points.empty() ? nullptr : points.front().data(),
                       static_cast<Py_ssize_t>(points.size()), 2);
}

bool readDoubles(PyObject *o, int cols, const char *what, std::vector<double> &out,
                 Py_ssize_t &rows) {
    Py_buffer view;
    PyObject *arr = contiguous(o, "float64", cols, what, view);
    if (!arr) return false;
    rows = view.shape[0];
    const double *p = static_cast<const double *>(view.buf);
    out.assign(p, p + view.len / static_cast<Py_ssize_t>(sizeof(double)));
    PyBuffer_Release(&view);
    Py_DECREF(arr);
    return true;
}

bool readIndices(PyObject *o, int cols, const char *what, std::vector<long long> &out,
                 Py_ssize_t &rows) {
    Py_buffer view;
    PyObject *arr = contiguous(o, "int64", cols, what, view);
    if (!arr) return false;
    rows = view.shape[0];
    const std::int64_t *p = static_cast<const std::int64_t *>(view.buf);
    out.assign(p, p + view.len / static_cast<Py_ssize_t>(sizeof(std::int64_t)));
    PyBuffer_Release(&view);
    Py_DECREF(arr);
    return true;
}

}  // namespace pycg
