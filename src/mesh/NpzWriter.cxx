#include "mesh/NpzWriter.hxx"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace {

// The descr numpy writes for each element type. The bytes go out in the host's
// order, which the '<' asserts is little-endian; hostIsLittleEndian() checks
// that once before anything is written.
template <class T> struct Descr;
template <> struct Descr<float> { static constexpr const char *value = "<f4"; };
template <> struct Descr<double> { static constexpr const char *value = "<f8"; };
template <> struct Descr<std::int8_t> { static constexpr const char *value = "|i1"; };
template <> struct Descr<std::int32_t> { static constexpr const char *value = "<i4"; };
template <> struct Descr<std::uint32_t> { static constexpr const char *value = "<u4"; };
template <> struct Descr<std::int64_t> { static constexpr const char *value = "<i8"; };

bool hostIsLittleEndian() {
    const std::uint16_t one = 1;
    unsigned char first;
    std::memcpy(&first, &one, 1);
    return first == 1;
}

// The .npy header of NEP 1, version 1.0: magic, version, a little-endian
// uint16 length, then a Python dict literal padded with spaces to end in a
// newline on a 64-byte boundary, which is what numpy itself writes so that
// the data that follows can be memory-mapped aligned.
std::vector<unsigned char> npyHeader(const std::string &descr, const std::vector<size_t> &shape) {
    std::string dims = "(";
    for (size_t i = 0; i < shape.size(); ++i) {
        dims += std::to_string(shape[i]);
        if (shape.size() == 1 || i + 1 < shape.size()) dims += ",";
        if (i + 1 < shape.size()) dims += " ";
    }
    dims += ")";
    std::string dict = "{'descr': '" + descr + "', 'fortran_order': False, 'shape': " + dims + ", }";
    const size_t fixed = 10;   // magic (6) + version (2) + length (2)
    size_t total = fixed + dict.size() + 1;
    total = (total + 63) / 64 * 64;
    dict.append(total - fixed - dict.size() - 1, ' ');
    dict.push_back('\n');
    if (dict.size() > std::numeric_limits<std::uint16_t>::max())
        throw std::invalid_argument("NpzWriter: .npy header too long");

    std::vector<unsigned char> h = {0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0};
    const std::uint16_t len = static_cast<std::uint16_t>(dict.size());
    h.push_back(static_cast<unsigned char>(len & 0xff));
    h.push_back(static_cast<unsigned char>(len >> 8));
    h.insert(h.end(), dict.begin(), dict.end());
    return h;
}

std::uint32_t crc32(const unsigned char *p, size_t n) {
    static std::uint32_t table[256];
    static bool built = false;
    if (!built) {
        for (std::uint32_t i = 0; i < 256; ++i) {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            table[i] = c;
        }
        built = true;
    }
    std::uint32_t c = 0xFFFFFFFFu;
    for (size_t i = 0; i < n; ++i) c = table[(c ^ p[i]) & 0xff] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

void put16(std::vector<unsigned char> &out, std::uint32_t v) {
    out.push_back(static_cast<unsigned char>(v & 0xff));
    out.push_back(static_cast<unsigned char>((v >> 8) & 0xff));
}

void put32(std::vector<unsigned char> &out, std::uint32_t v) {
    put16(out, v & 0xffff);
    put16(out, v >> 16);
}

// UTF-8 to code points, for numpy's <U (UTF-32). A malformed byte becomes
// U+FFFD rather than an error: the strings written are names and labels.
std::vector<std::uint32_t> codePoints(const std::string &s) {
    std::vector<std::uint32_t> out;
    for (size_t i = 0; i < s.size();) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        const int extra = c < 0x80 ? 0 : (c >> 5) == 0x6 ? 1 : (c >> 4) == 0xE ? 2 : (c >> 3) == 0x1E ? 3 : -1;
        if (extra < 0 || i + extra >= s.size() + (extra == 0)) {
            out.push_back(0xFFFD);
            ++i;
            continue;
        }
        std::uint32_t cp = extra == 0 ? c : c & (0x3F >> extra);
        bool ok = true;
        for (int k = 1; k <= extra; ++k) {
            const unsigned char d = static_cast<unsigned char>(s[i + k]);
            if ((d >> 6) != 0x2) { ok = false; break; }
            cp = (cp << 6) | (d & 0x3F);
        }
        out.push_back(ok ? cp : 0xFFFD);
        i += ok ? extra + 1 : 1;
    }
    return out;
}

}  // namespace

void NpzWriter::addRaw(const std::string &name, const std::string &descr, const std::vector<size_t> &shape,
                       const void *data, size_t bytes) {
    Member m;
    m.name = name;
    m.bytes = npyHeader(descr, shape);
    const unsigned char *p = static_cast<const unsigned char *>(data);
    m.bytes.insert(m.bytes.end(), p, p + bytes);
    for (Member &old : members_) {
        if (old.name == name) {
            old = std::move(m);
            return;
        }
    }
    members_.push_back(std::move(m));
}

template <class T>
void NpzWriter::add(const std::string &name, const std::vector<size_t> &shape, const std::vector<T> &data) {
    size_t count = 1;
    for (size_t d : shape) count *= d;
    if (count != data.size())
        throw std::invalid_argument("NpzWriter: '" + name + "' has " + std::to_string(data.size()) +
                                    " values for a shape of " + std::to_string(count));
    addRaw(name, Descr<T>::value, shape, data.data(), data.size() * sizeof(T));
}

template void NpzWriter::add<float>(const std::string &, const std::vector<size_t> &, const std::vector<float> &);
template void NpzWriter::add<double>(const std::string &, const std::vector<size_t> &, const std::vector<double> &);
template void NpzWriter::add<std::int8_t>(const std::string &, const std::vector<size_t> &,
                                          const std::vector<std::int8_t> &);
template void NpzWriter::add<std::int32_t>(const std::string &, const std::vector<size_t> &,
                                           const std::vector<std::int32_t> &);
template void NpzWriter::add<std::uint32_t>(const std::string &, const std::vector<size_t> &,
                                            const std::vector<std::uint32_t> &);
template void NpzWriter::add<std::int64_t>(const std::string &, const std::vector<size_t> &,
                                           const std::vector<std::int64_t> &);

void NpzWriter::addScalar(const std::string &name, double value) {
    add<double>(name, {}, std::vector<double>{value});
}

void NpzWriter::addScalar(const std::string &name, long long value) {
    add<std::int64_t>(name, {}, std::vector<std::int64_t>{static_cast<std::int64_t>(value)});
}

void NpzWriter::addString(const std::string &name, const std::string &utf8) {
    std::vector<std::uint32_t> cp = codePoints(utf8);
    // numpy has no <U0; an empty string is stored as <U1 holding a NUL, which
    // numpy reads back as ''.
    if (cp.empty()) cp.push_back(0);
    addRaw(name, "<U" + std::to_string(cp.size()), {}, cp.data(), cp.size() * sizeof(std::uint32_t));
}

bool NpzWriter::has(const std::string &name) const {
    return std::any_of(members_.begin(), members_.end(), [&](const Member &m) { return m.name == name; });
}

// The zip container (APPNOTE 6.3.x, the subset every reader accepts): per
// member a local file header and its bytes, then the central directory and
// its end record. Method 0 (stored), no data descriptors, no Zip64.
bool NpzWriter::write(const std::string &path) const {
    if (!hostIsLittleEndian()) return false;
    const std::uint32_t kMax = std::numeric_limits<std::uint32_t>::max();
    const std::uint32_t dosTime = 0, dosDate = (0u << 9) | (1u << 5) | 1u;   // 1980-01-01 00:00

    std::vector<unsigned char> central;
    std::vector<std::uint32_t> offsets;
    size_t offset = 0;

    FILE *f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    bool ok = true;
    for (const Member &m : members_) {
        const std::string file = m.name + ".npy";
        if (m.bytes.size() >= kMax || offset >= kMax || file.size() > 0xffff) { ok = false; break; }
        const std::uint32_t crc = crc32(m.bytes.data(), m.bytes.size());
        const std::uint32_t size = static_cast<std::uint32_t>(m.bytes.size());

        std::vector<unsigned char> local;
        put32(local, 0x04034b50);
        put16(local, 20);            // version needed: 2.0
        put16(local, 0);             // flags
        put16(local, 0);             // stored
        put16(local, dosTime);
        put16(local, dosDate);
        put32(local, crc);
        put32(local, size);
        put32(local, size);
        put16(local, static_cast<std::uint32_t>(file.size()));
        put16(local, 0);             // extra field length
        local.insert(local.end(), file.begin(), file.end());

        put32(central, 0x02014b50);
        put16(central, 20);          // version made by
        put16(central, 20);
        put16(central, 0);
        put16(central, 0);
        put16(central, dosTime);
        put16(central, dosDate);
        put32(central, crc);
        put32(central, size);
        put32(central, size);
        put16(central, static_cast<std::uint32_t>(file.size()));
        put16(central, 0);           // extra
        put16(central, 0);           // comment
        put16(central, 0);           // disk
        put16(central, 0);           // internal attributes
        put32(central, 0);           // external attributes
        put32(central, static_cast<std::uint32_t>(offset));
        central.insert(central.end(), file.begin(), file.end());

        ok = std::fwrite(local.data(), 1, local.size(), f) == local.size() &&
             std::fwrite(m.bytes.data(), 1, m.bytes.size(), f) == m.bytes.size();
        if (!ok) break;
        offset += local.size() + m.bytes.size();
    }
    if (ok && (offset >= kMax || members_.size() > 0xffff)) ok = false;
    if (ok) {
        std::vector<unsigned char> end;
        put32(end, 0x06054b50);
        put16(end, 0);               // this disk
        put16(end, 0);               // disk with the directory
        put16(end, static_cast<std::uint32_t>(members_.size()));
        put16(end, static_cast<std::uint32_t>(members_.size()));
        put32(end, static_cast<std::uint32_t>(central.size()));
        put32(end, static_cast<std::uint32_t>(offset));
        put16(end, 0);               // comment length
        ok = std::fwrite(central.data(), 1, central.size(), f) == central.size() &&
             std::fwrite(end.data(), 1, end.size(), f) == end.size();
    }
    if (std::fclose(f) != 0) ok = false;
    if (!ok) std::remove(path.c_str());
    return ok;
}
