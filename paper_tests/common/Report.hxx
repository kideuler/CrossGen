#ifndef __PAPER_REPORT_HXX__
#define __PAPER_REPORT_HXX__

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

// The console and the CSV, so that no experiment in this directory has to grow
// its own printing. Every experiment writes both: the console line is what you
// read while it runs, and the CSV is what the plots and the per-model
// supplement are built from.
namespace paper {

inline void ensureDir(const std::string &path) {
    std::string cur;
    for (std::size_t i = 0; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            if (!cur.empty()) ::mkdir(cur.c_str(), 0755);
        }
        if (i < path.size()) cur.push_back(path[i]);
    }
}

// A CSV with a fixed header. Rows are written as they are produced, so a run
// that is interrupted still leaves everything it had finished.
class Csv {
public:
    Csv(const std::string &path, std::vector<std::string> header)
        : out(path), cols(std::move(header)) {
        if (!out) {
            std::cerr << "paper_tests: could not open " << path << " for writing\n";
            return;
        }
        for (std::size_t i = 0; i < cols.size(); ++i) out << (i ? "," : "") << cols[i];
        out << "\n";
        out << std::setprecision(10);
    }

    void row(const std::vector<std::string> &values) {
        for (std::size_t i = 0; i < values.size(); ++i) out << (i ? "," : "") << values[i];
        out << "\n";
        out.flush();
    }

    // Named form, so that a caller adding a column in one place cannot silently
    // shift every value in the row.
    void row(const std::map<std::string, std::string> &values) {
        std::vector<std::string> v;
        v.reserve(cols.size());
        for (const std::string &c : cols) {
            auto it = values.find(c);
            v.push_back(it == values.end() ? "" : it->second);
        }
        row(v);
    }

    bool good() const { return static_cast<bool>(out); }

private:
    std::ofstream out;
    std::vector<std::string> cols;
};

inline std::string num(double x, int digits = 6) {
    std::ostringstream oss;
    oss << std::setprecision(digits) << x;
    return oss.str();
}
inline std::string num(int x) { return std::to_string(x); }
inline std::string num(bool x) { return x ? "1" : "0"; }

// --- console ---------------------------------------------------------------

inline const char *kPass = "\033[32m[PASS]\033[0m";
inline const char *kFail = "\033[31m[FAIL]\033[0m";
inline const char *kWarn = "\033[33m[WARN]\033[0m";

struct Verdicts {
    int failures = 0;
    int warnings = 0;

    void check(bool ok, const std::string &what) {
        std::cout << "  " << (ok ? kPass : kFail) << " " << what << "\n";
        if (!ok) ++failures;
    }
    void warn(const std::string &what) {
        std::cout << "  " << kWarn << " " << what << "\n";
        ++warnings;
    }
    void note(const std::string &what) { std::cout << "         " << what << "\n"; }
};

inline void heading(const std::string &title) {
    std::cout << "\n" << title << "\n" << std::string(title.size(), '-') << "\n";
}

inline void banner(const std::string &title) {
    std::cout << "\n" << std::string(78, '=') << "\n" << title << "\n"
              << std::string(78, '=') << "\n";
}

// A fixed-width console table. Columns size themselves to the widest cell.
class Table {
public:
    explicit Table(std::vector<std::string> header) : rows{std::move(header)} {}
    void row(std::vector<std::string> r) { rows.push_back(std::move(r)); }
    void print(std::ostream &os = std::cout, const std::string &indent = "  ") const {
        if (rows.empty()) return;
        std::size_t n = 0;
        for (const auto &r : rows) n = std::max(n, r.size());
        std::vector<std::size_t> w(n, 0);
        for (const auto &r : rows)
            for (std::size_t i = 0; i < r.size(); ++i) w[i] = std::max(w[i], r[i].size());
        for (std::size_t k = 0; k < rows.size(); ++k) {
            os << indent;
            for (std::size_t i = 0; i < n; ++i) {
                const std::string &c = i < rows[k].size() ? rows[k][i] : std::string();
                os << (i == 0 ? "" : "  ") << std::left << std::setw(static_cast<int>(w[i])) << c;
            }
            os << "\n";
            if (k == 0) {
                os << indent;
                for (std::size_t i = 0; i < n; ++i)
                    os << (i == 0 ? "" : "  ") << std::string(w[i], '-');
                os << "\n";
            }
        }
        os << std::right;
    }

private:
    std::vector<std::vector<std::string>> rows;
};

} // namespace paper

#endif // __PAPER_REPORT_HXX__
