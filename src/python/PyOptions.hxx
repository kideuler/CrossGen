#ifndef __PY_OPTIONS_HXX__
#define __PY_OPTIONS_HXX__

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <initializer_list>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// Keyword arguments onto a C++ Options struct, and a Report struct back out as
// a dict, for the crossgen module (src/python/CrossGenModule.cxx).
//
// Every method in this codebase is already controlled by one plain struct --
// ZIPLINE::Options, MERIDIAN::Options, mesh::TMOP::Options and so on -- whose
// defaults are the ones every corpus number was measured at. The binding does
// not restate those defaults, and it does not pick a subset of the fields to
// expose. A Table names each field once, by the path the C++ spells it by, and
// the Python keyword is derived from that path:
//
//     CG_OPTION(ZIPLINE::Options, trace.maxStepsPerSeparatrix)
//         -> keyword  trace_max_steps_per_separatrix
//
// so m.zipline(trace_max_steps_per_separatrix=5000) sets exactly that member,
// and anything not passed keeps the value the struct's own initialiser gave it.
// The name is computed from the member path rather than typed a second time,
// which is what keeps the two from drifting apart when a field is renamed: the
// table then stops compiling instead of silently answering to the old name.
//
// A field of a type with no conversion below is a compile error in the table,
// not a silent omission, so a table either covers a member faithfully or not
// at all.
namespace pycg {

// camelCase member path -> snake_case keyword: `.` becomes `_`, and an
// underscore goes in front of every upper-case letter that starts a word
// ("dualMBOGamma" -> "dual_mbo_gamma", "halfOGrids" -> "half_o_grids",
// "q5Verified" -> "q5_verified").
std::string snakeCase(const char *path);

// ---------------------------------------------------------------------------
// Enumerations
//
// EnumNames<E>::list is E's {value, keyword} table, names in snake_case. It is
// specialised next to the table that first needs it (PyTables.hxx for the ones
// several methods share), so this header depends on no pipeline. A Python
// caller may give either the name or the integer value; a dict always shows
// the name.
// ---------------------------------------------------------------------------
template <class E>
struct EnumNames;

// A TypeError naming the valid choices. Always returns -1.
int enumError(PyObject *value, const std::vector<const char *> &names);

// ---------------------------------------------------------------------------
// Conversions. toPy returns a new reference, or null with an exception set.
// fromPy returns 0, or -1 with a TypeError/ValueError set.
//
// They are strict where the C++ is: an int field refuses 2.5 rather than
// truncating it, and a bool field takes True/False (or 0/1) and nothing else,
// because a caller who passed a float where a count was meant has made a
// mistake the struct would otherwise hide.
// ---------------------------------------------------------------------------
PyObject *toPy(bool v);
PyObject *toPy(int v);
PyObject *toPy(unsigned v);
PyObject *toPy(long long v);
PyObject *toPy(unsigned long long v);
PyObject *toPy(double v);
PyObject *toPy(const std::string &v);
PyObject *toPy(const std::vector<double> &v);
PyObject *toPy(const std::vector<std::string> &v);

int fromPy(PyObject *o, bool &v);
int fromPy(PyObject *o, int &v);
int fromPy(PyObject *o, unsigned &v);
int fromPy(PyObject *o, unsigned long long &v);
int fromPy(PyObject *o, double &v);
int fromPy(PyObject *o, std::vector<double> &v);

template <class E, std::enable_if_t<std::is_enum<E>::value, int> = 0>
PyObject *toPy(E v) {
    for (const auto &entry : EnumNames<E>::list)
        if (entry.first == v) return PyUnicode_FromString(entry.second);
    // A value the table has no name for: say what it is rather than fail.
    return PyLong_FromLongLong(static_cast<long long>(v));
}

template <class E, std::enable_if_t<std::is_enum<E>::value, int> = 0>
int fromPy(PyObject *o, E &v) {
    std::vector<const char *> names;
    for (const auto &entry : EnumNames<E>::list) names.push_back(entry.second);
    if (PyUnicode_Check(o)) {
        const char *s = PyUnicode_AsUTF8(o);
        if (!s) return -1;
        for (const auto &entry : EnumNames<E>::list) {
            if (std::string(s) == entry.second) {
                v = entry.first;
                return 0;
            }
        }
        return enumError(o, names);
    }
    if (PyLong_Check(o) && !PyBool_Check(o)) {
        const long long k = PyLong_AsLongLong(o);
        if (k == -1 && PyErr_Occurred()) return -1;
        for (const auto &entry : EnumNames<E>::list) {
            if (static_cast<long long>(entry.first) == k) {
                v = entry.first;
                return 0;
            }
        }
    }
    return enumError(o, names);
}

// ---------------------------------------------------------------------------
// Fields and tables
// ---------------------------------------------------------------------------

// One member of O: its keyword, a getter, and a setter (null for a report
// field, which Python may read and not write).
template <class O>
struct Field {
    std::string name;
    PyObject *(*get)(const O &);
    int (*set)(O &, PyObject *);

    Field(const char *path, PyObject *(*g)(const O &), int (*s)(O &, PyObject *))
        : name(snakeCase(path)), get(g), set(s) {}
};

// The type-erased half, so that one call can spread a single **kwargs over
// several structs of different types (BlockDecomposition.mesh takes the
// mesher's options and mesh::QuadMesh's node options in the same call).
class TableBase {
public:
    virtual ~TableBase() = default;

    const std::vector<std::string> &names() const { return names_; }
    // Index of the field called `name`, or -1.
    int find(const std::string &name) const;

    virtual PyObject *get(const void *object, int field) const = 0;
    // -1 with an exception set on a bad value, or on a read-only field.
    virtual int set(void *object, int field, PyObject *value) const = 0;

protected:
    // Throws std::logic_error if two fields snake-case to one keyword: a table
    // that cannot tell two members apart is a defect in the table, and the
    // module refuses to import rather than set the wrong one.
    void index(std::vector<std::string> names);

private:
    std::vector<std::string> names_;
    std::unordered_map<std::string, int> byName_;
};

template <class O>
class Table final : public TableBase {
public:
    Table(std::initializer_list<Field<O>> fields) : fields_(fields) {
        std::vector<std::string> names;
        names.reserve(fields_.size());
        for (const Field<O> &f : fields_) names.push_back(f.name);
        index(std::move(names));
    }

    PyObject *get(const void *object, int field) const override {
        return fields_[field].get(*static_cast<const O *>(object));
    }

    int set(void *object, int field, PyObject *value) const override {
        const Field<O> &f = fields_[field];
        if (!f.set) {
            PyErr_Format(PyExc_AttributeError, "'%s' is read-only", f.name.c_str());
            return -1;
        }
        return f.set(*static_cast<O *>(object), value);
    }

private:
    std::vector<Field<O>> fields_;
};

// A struct to write keyword arguments into, and a struct to read a dict out
// of. `prefix` goes in front of every key of that struct in the dict, for a
// report assembled from several stages' reports whose field names overlap.
struct Target {
    const TableBase *table;
    void *object;
};
struct Source {
    const TableBase *table;
    const void *object;
    const char *prefix;
};

template <class O>
Target target(const Table<O> &t, O &o) { return Target{&t, &o}; }
template <class O>
Source source(const Table<O> &t, const O &o, const char *prefix = "") {
    return Source{&t, &o, prefix};
}

// Sets each keyword argument of `kwargs` (a dict, or null for none) on the
// first target whose table has a field of that name. An unknown keyword is a
// TypeError that names `function` and suggests the nearest valid keywords.
// Returns 0, or -1 with an exception set.
int applyKwargs(PyObject *kwargs, const char *function, std::initializer_list<Target> targets);

// A dict of every field of every source, in table order. New reference.
PyObject *toDict(std::initializer_list<Source> sources);
// The same into an existing dict. Returns 0, or -1 with an exception set.
int addToDict(PyObject *dict, const Source &src);

}  // namespace pycg

// One member of an Options struct: read and written from Python.
#define CG_OPTION(O, path)                                                      \
    ::pycg::Field<O>(                                                           \
        #path, [](const O &o) -> PyObject * { return ::pycg::toPy(o.path); },   \
        [](O &o, PyObject *v) -> int { return ::pycg::fromPy(v, o.path); })

// One member of a Report or Status struct: read from Python only.
#define CG_REPORT(O, path)                                                      \
    ::pycg::Field<O>(                                                           \
        #path, [](const O &o) -> PyObject * { return ::pycg::toPy(o.path); },   \
        nullptr)

#endif // __PY_OPTIONS_HXX__
