#ifndef __PY_UTIL_HXX__
#define __PY_UTIL_HXX__

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <array>
#include <exception>
#include <mutex>
#include <new>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// Arrays in and out of the crossgen module, and the one rule for running C++
// work from it.
//
// ### Arrays
//
// Everything array-shaped comes back as a numpy.ndarray -- float64 for
// coordinates and eigenvalues, int64 for indices and material ids -- and is a
// copy: writing into it changes nothing on the C++ side. numpy is imported the
// first time an array is made rather than linked against at build time, so the
// module builds against Python's headers alone and works with whichever numpy
// the interpreter has.
//
// ### Running C++ work
//
// Every pipeline, mesher, smoother and eigensolver call goes through
// compute(), which does three things in a fixed order:
//
//   1. releases the GIL, so other Python threads run while a pipeline does;
//   2. takes one process-wide lock, so no two crossgen computations run at
//      once -- Triangle, which ATLAS's coarse domains and others mesh through,
//      is not re-entrant (see ATLAS::run), and the pipelines were never written
//      to share a process with themselves. The parallelism is inside each call,
//      on OpenMP and ATLAS's own search threads; for many models at once, use
//      processes;
//   3. turns any C++ exception into a Python RuntimeError (MemoryError for
//      std::bad_alloc) once the GIL is back.
//
// The GIL is always given up before the lock is taken and taken back after the
// lock is released, which is what keeps a thread waiting on one from ever
// holding the other.
namespace pycg {

// ---- numpy out --------------------------------------------------------------
// New references, or null with an exception set.
PyObject *doubleArray(const double *data, Py_ssize_t rows, Py_ssize_t cols);   // cols 0: 1-D
PyObject *indexArray(const int *data, Py_ssize_t rows, Py_ssize_t cols);       // int64
PyObject *doubleArray(const std::vector<double> &v);
PyObject *indexArray(const std::vector<int> &v);
PyObject *pointArray(const std::vector<Point> &points);                        // (n, 2)

template <size_t K>
PyObject *indexArray(const std::vector<std::array<int, K>> &rows) {
    return indexArray(rows.empty() ? nullptr : rows.front().data(),
                      static_cast<Py_ssize_t>(rows.size()), static_cast<Py_ssize_t>(K));
}

// ---- numpy in ---------------------------------------------------------------
// `o` read through numpy.ascontiguousarray as float64 of shape (n, cols) -- or
// (n,) when cols is 0 -- into `out`, row-major. `what` names the argument in
// the error. Return false with an exception set.
bool readDoubles(PyObject *o, int cols, const char *what, std::vector<double> &out, Py_ssize_t &rows);
bool readIndices(PyObject *o, int cols, const char *what, std::vector<long long> &out, Py_ssize_t &rows);

// ---- running C++ ------------------------------------------------------------
std::mutex &computeMutex();

class ReleasedGIL {
public:
    ReleasedGIL() : state_(PyEval_SaveThread()) {}
    ~ReleasedGIL() { PyEval_RestoreThread(state_); }
    ReleasedGIL(const ReleasedGIL &) = delete;
    ReleasedGIL &operator=(const ReleasedGIL &) = delete;

private:
    PyThreadState *state_;
};

// Run `f` as described above. `f` must not touch any Python object. Returns
// false with a Python exception set if `f` threw.
template <class F>
bool compute(F &&f) {
    std::string error;
    bool failed = false, noMemory = false;
    {
        ReleasedGIL released;
        std::lock_guard<std::mutex> lock(computeMutex());
        try {
            f();
        } catch (const std::bad_alloc &) {
            failed = noMemory = true;
        } catch (const std::exception &e) {
            failed = true;
            error = e.what();
        } catch (...) {
            failed = true;
            error = "unknown C++ exception";
        }
    }
    if (!failed) return true;
    if (noMemory) PyErr_NoMemory();
    else PyErr_SetString(PyExc_RuntimeError, error.c_str());
    return false;
}

}  // namespace pycg

#endif // __PY_UTIL_HXX__
