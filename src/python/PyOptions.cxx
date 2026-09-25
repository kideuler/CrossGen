#include "python/PyOptions.hxx"

#include <cctype>
#include <climits>
#include <stdexcept>

namespace pycg {

std::string snakeCase(const char *path) {
    std::string out;
    const std::string s(path);
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c == '.') {
            out += '_';
            continue;
        }
        if (std::isupper(c)) {
            // A word starts here if the letter before it ends one ("dualMBO":
            // the M) or if this is the last capital of an acronym running into
            // a lower-case word ("MBOGamma": the G). At the start of the path,
            // or straight after a dot, it is already a word boundary.
            const bool first = (i == 0 || s[i - 1] == '.');
            const unsigned char prev = first ? 0 : static_cast<unsigned char>(s[i - 1]);
            const unsigned char next = (i + 1 < s.size()) ? static_cast<unsigned char>(s[i + 1]) : 0;
            const bool boundary =
                !first && (std::islower(prev) || std::isdigit(prev) ||
                           (std::isupper(prev) && std::islower(next)));
            if (boundary) out += '_';
            out += static_cast<char>(std::tolower(c));
            continue;
        }
        out += static_cast<char>(c);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Conversions
// ---------------------------------------------------------------------------
PyObject *toPy(bool v) { return PyBool_FromLong(v ? 1 : 0); }
PyObject *toPy(int v) { return PyLong_FromLong(v); }
PyObject *toPy(unsigned v) { return PyLong_FromUnsignedLong(v); }
PyObject *toPy(long long v) { return PyLong_FromLongLong(v); }
PyObject *toPy(unsigned long long v) { return PyLong_FromUnsignedLongLong(v); }
PyObject *toPy(double v) { return PyFloat_FromDouble(v); }

PyObject *toPy(const std::string &v) {
    // The pipelines' messages are UTF-8 prose; "replace" keeps a stray byte
    // from turning a report into an exception.
    return PyUnicode_DecodeUTF8(v.data(), static_cast<Py_ssize_t>(v.size()), "replace");
}

PyObject *toPy(const std::vector<double> &v) {
    PyObject *list = PyList_New(static_cast<Py_ssize_t>(v.size()));
    if (!list) return nullptr;
    for (size_t i = 0; i < v.size(); ++i) {
        PyObject *x = PyFloat_FromDouble(v[i]);
        if (!x) { Py_DECREF(list); return nullptr; }
        PyList_SET_ITEM(list, static_cast<Py_ssize_t>(i), x);
    }
    return list;
}

PyObject *toPy(const std::vector<std::string> &v) {
    PyObject *list = PyList_New(static_cast<Py_ssize_t>(v.size()));
    if (!list) return nullptr;
    for (size_t i = 0; i < v.size(); ++i) {
        PyObject *x = toPy(v[i]);
        if (!x) { Py_DECREF(list); return nullptr; }
        PyList_SET_ITEM(list, static_cast<Py_ssize_t>(i), x);
    }
    return list;
}

namespace {

// A Python integer, or anything with __index__ (numpy's integer scalars are not
// PyLong subclasses but do have it). A bool is refused: it is an int to Python
// and almost never what a count was meant to be.
int readIndex(PyObject *o, long long lo, long long hi, long long &out) {
    if (PyBool_Check(o) || !PyIndex_Check(o)) {
        PyErr_Format(PyExc_TypeError, "expected an integer, got %.100s", Py_TYPE(o)->tp_name);
        return -1;
    }
    PyObject *i = PyNumber_Index(o);
    if (!i) return -1;
    int overflow = 0;
    const long long k = PyLong_AsLongLongAndOverflow(i, &overflow);
    Py_DECREF(i);
    if (k == -1 && PyErr_Occurred()) return -1;
    if (overflow != 0 || k < lo || k > hi) {
        PyErr_Format(PyExc_OverflowError, "integer %R out of range [%lld, %lld]", o, lo, hi);
        return -1;
    }
    out = k;
    return 0;
}

}  // namespace

int fromPy(PyObject *o, bool &v) {
    if (PyBool_Check(o)) {
        v = (o == Py_True);
        return 0;
    }
    if (PyLong_Check(o)) {
        const long k = PyLong_AsLong(o);
        if (k == 0 || k == 1) {
            v = (k == 1);
            return 0;
        }
    }
    PyErr_Format(PyExc_TypeError, "expected True or False, got %R", o);
    return -1;
}

int fromPy(PyObject *o, int &v) {
    long long k = 0;
    if (readIndex(o, INT_MIN, INT_MAX, k) < 0) return -1;
    v = static_cast<int>(k);
    return 0;
}

int fromPy(PyObject *o, unsigned &v) {
    long long k = 0;
    if (readIndex(o, 0, UINT_MAX, k) < 0) return -1;
    v = static_cast<unsigned>(k);
    return 0;
}

int fromPy(PyObject *o, unsigned long long &v) {
    if (PyBool_Check(o) || !PyIndex_Check(o)) {
        PyErr_Format(PyExc_TypeError, "expected an integer, got %.100s", Py_TYPE(o)->tp_name);
        return -1;
    }
    PyObject *i = PyNumber_Index(o);
    if (!i) return -1;
    const unsigned long long k = PyLong_AsUnsignedLongLong(i);
    Py_DECREF(i);
    if (PyErr_Occurred()) return -1;
    v = k;
    return 0;
}

int fromPy(PyObject *o, double &v) {
    if (PyBool_Check(o)) {
        PyErr_SetString(PyExc_TypeError, "expected a number, got a bool");
        return -1;
    }
    const double x = PyFloat_AsDouble(o);
    if (x == -1.0 && PyErr_Occurred()) return -1;
    v = x;
    return 0;
}

int fromPy(PyObject *o, std::vector<double> &v) {
    if (PyUnicode_Check(o) || PyBytes_Check(o)) {
        PyErr_SetString(PyExc_TypeError, "expected a sequence of numbers, got a string");
        return -1;
    }
    PyObject *seq = PySequence_Fast(o, "expected a sequence of numbers");
    if (!seq) return -1;
    const Py_ssize_t n = PySequence_Fast_GET_SIZE(seq);
    std::vector<double> out(static_cast<size_t>(n));
    for (Py_ssize_t i = 0; i < n; ++i) {
        if (fromPy(PySequence_Fast_GET_ITEM(seq, i), out[static_cast<size_t>(i)]) < 0) {
            Py_DECREF(seq);
            return -1;
        }
    }
    Py_DECREF(seq);
    v.swap(out);
    return 0;
}

int enumError(PyObject *value, const std::vector<const char *> &names) {
    std::string list;
    for (size_t i = 0; i < names.size(); ++i) {
        if (i) list += ", ";
        list += "'";
        list += names[i];
        list += "'";
    }
    PyObject *type = (PyUnicode_Check(value) || PyLong_Check(value)) ? PyExc_ValueError
                                                                     : PyExc_TypeError;
    PyErr_Format(type, "expected one of %s (or the integer value), got %R", list.c_str(), value);
    return -1;
}

// ---------------------------------------------------------------------------
// Tables
// ---------------------------------------------------------------------------
void TableBase::index(std::vector<std::string> names) {
    names_ = std::move(names);
    byName_.clear();
    for (size_t i = 0; i < names_.size(); ++i) {
        if (!byName_.emplace(names_[i], static_cast<int>(i)).second)
            throw std::logic_error("two fields map to the keyword '" + names_[i] + "'");
    }
}

int TableBase::find(const std::string &name) const {
    const auto it = byName_.find(name);
    return it == byName_.end() ? -1 : it->second;
}

namespace {

// Put `prefix` in front of the message of the exception being raised, keeping
// its type: "expected an integer" becomes "zipline(): trace_max_steps_per_
// separatrix: expected an integer".
void prefixError(const std::string &prefix) {
#if PY_VERSION_HEX >= 0x030C0000
    PyObject *value = PyErr_GetRaisedException();
    PyObject *type = value ? reinterpret_cast<PyObject *>(Py_TYPE(value)) : PyExc_TypeError;
    Py_INCREF(type);
#else
    PyObject *type = nullptr, *value = nullptr, *trace = nullptr;
    PyErr_Fetch(&type, &value, &trace);
    PyErr_NormalizeException(&type, &value, &trace);
    Py_XDECREF(trace);
    if (!type) { type = PyExc_TypeError; Py_INCREF(type); }
#endif
    PyObject *text = value ? PyObject_Str(value) : nullptr;
    const char *msg = text ? PyUnicode_AsUTF8(text) : nullptr;
    PyErr_Clear();
    PyErr_Format(type, "%s: %s", prefix.c_str(), msg ? msg : "invalid value");
    Py_XDECREF(text);
    Py_DECREF(type);
    Py_XDECREF(value);
}

// difflib's nearest matches to `key` among `names`, as "'a', 'b' or 'c'", or
// empty. A failure here is swallowed: a missing suggestion is no reason to lose
// the real error.
std::string suggestions(const std::string &key, const std::vector<std::string> &names) {
    std::string out;
    PyObject *difflib = PyImport_ImportModule("difflib");
    PyObject *list = PyList_New(0);
    PyObject *matches = nullptr;
    if (difflib && list) {
        for (const std::string &n : names) {
            PyObject *s = PyUnicode_FromString(n.c_str());
            if (s) { PyList_Append(list, s); Py_DECREF(s); }
        }
        matches = PyObject_CallMethod(difflib, "get_close_matches", "sOid", key.c_str(), list, 3, 0.6);
    }
    if (matches && PyList_Check(matches)) {
        const Py_ssize_t n = PyList_GET_SIZE(matches);
        for (Py_ssize_t i = 0; i < n; ++i) {
            const char *m = PyUnicode_AsUTF8(PyList_GET_ITEM(matches, i));
            if (!m) continue;
            if (i > 0) out += (i + 1 == n) ? " or " : ", ";
            out += "'";
            out += m;
            out += "'";
        }
    }
    PyErr_Clear();
    Py_XDECREF(matches);
    Py_XDECREF(list);
    Py_XDECREF(difflib);
    return out;
}

}  // namespace

int applyKwargs(PyObject *kwargs, const char *function, std::initializer_list<Target> targets) {
    if (!kwargs) return 0;
    PyObject *key = nullptr, *value = nullptr;
    Py_ssize_t pos = 0;
    while (PyDict_Next(kwargs, &pos, &key, &value)) {
        const char *k = PyUnicode_AsUTF8(key);
        if (!k) return -1;
        bool found = false;
        for (const Target &t : targets) {
            const int f = t.table->find(k);
            if (f < 0) continue;
            found = true;
            if (t.table->set(t.object, f, value) < 0) {
                prefixError(std::string(function) + "(): " + k);
                return -1;
            }
            break;
        }
        if (!found) {
            std::vector<std::string> all;
            for (const Target &t : targets)
                all.insert(all.end(), t.table->names().begin(), t.table->names().end());
            const std::string near = suggestions(k, all);
            PyErr_Format(PyExc_TypeError,
                         "%s() got an unexpected keyword argument '%s'%s%s; crossgen.options() "
                         "lists the valid ones",
                         function, k, near.empty() ? "" : " -- did you mean ", near.c_str());
            return -1;
        }
    }
    return 0;
}

int addToDict(PyObject *dict, const Source &src) {
    const std::vector<std::string> &names = src.table->names();
    for (size_t i = 0; i < names.size(); ++i) {
        PyObject *v = src.table->get(src.object, static_cast<int>(i));
        if (!v) return -1;
        const std::string key = std::string(src.prefix ? src.prefix : "") + names[i];
        const int rc = PyDict_SetItemString(dict, key.c_str(), v);
        Py_DECREF(v);
        if (rc < 0) return -1;
    }
    return 0;
}

PyObject *toDict(std::initializer_list<Source> sources) {
    PyObject *dict = PyDict_New();
    if (!dict) return nullptr;
    for (const Source &s : sources) {
        if (addToDict(dict, s) < 0) {
            Py_DECREF(dict);
            return nullptr;
        }
    }
    return dict;
}

}  // namespace pycg
