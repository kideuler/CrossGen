#ifndef __NPZ_WRITER_HXX__
#define __NPZ_WRITER_HXX__

#include <cstddef>
#include <string>
#include <vector>

// A numpy .npz archive written from C++: what numpy.savez writes and
// numpy.load reads, so that data made here reaches Python without a second
// format in between and without a converter to keep in step with it.
//
// An .npz is a zip archive whose members are .npy files, one per array, each
// a short text header (dtype, shape, order) followed by the raw C-order bytes
// (numpy's NEP 1, format version 1.0). The members are *stored*, not deflated:
// numpy.load accepts either, stored needs no zlib, and the arrays this writes
// -- coordinates, barycentric weights, indices -- are float and integer data
// that deflate would shrink by little while costing every training-time load
// the time to inflate them. No member, and no archive, may exceed 4 GiB (no
// Zip64); write() refuses rather than produce an archive numpy cannot read.
//
// The array keys are the member names without ".npy", as numpy.savez would
// have them. Adding a key twice keeps the second.
class NpzWriter {
public:
    // An array of shape `shape` (empty: a 0-d scalar) from `data`, C order.
    // T is float, double, int8_t, int32_t, uint32_t or int64_t, which become
    // numpy's <f4, <f8, |i1, <i4, <u4 and <i8. Throws std::invalid_argument if
    // the shape's product is not data.size().
    template <class T>
    void add(const std::string &name, const std::vector<size_t> &shape, const std::vector<T> &data);

    // 0-d scalars and a 0-d unicode string (numpy <U, as numpy.str_ reads back).
    void addScalar(const std::string &name, double value);
    void addScalar(const std::string &name, long long value);
    void addString(const std::string &name, const std::string &utf8);

    bool has(const std::string &name) const;

    // False if the file cannot be written or a member is too large for a
    // plain zip.
    bool write(const std::string &path) const;

private:
    struct Member {
        std::string name;              // without ".npy"
        std::vector<unsigned char> bytes;
    };
    void addRaw(const std::string &name, const std::string &descr, const std::vector<size_t> &shape,
                const void *data, size_t bytes);

    std::vector<Member> members_;
};

#endif // __NPZ_WRITER_HXX__
