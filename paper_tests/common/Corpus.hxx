#ifndef __PAPER_CORPUS_HXX__
#define __PAPER_CORPUS_HXX__

#include <algorithm>
#include <dirent.h>
#include <memory>
#include <string>
#include <vector>

#include "mesh/Mesh.hxx"

// Where data/meshes is. Set by CMake so an experiment can be run from anywhere.
#ifndef PAPER_MESH_DIR
#define PAPER_MESH_DIR "data/meshes"
#endif

namespace paper {

struct Model {
    std::string name;   // "geom007"
    std::string path;
    std::shared_ptr<Mesh> mesh;   // loaded lazily by `load`
};

// Every .obj in data/meshes/<sub>, in name order.
inline std::vector<Model> corpus(const std::string &sub,
                                 const std::string &root = PAPER_MESH_DIR) {
    std::vector<Model> out;
    const std::string dir = root + "/" + sub;
    DIR *d = ::opendir(dir.c_str());
    if (!d) return out;
    while (dirent *e = ::readdir(d)) {
        const std::string n = e->d_name;
        if (n.size() < 5 || n.substr(n.size() - 4) != ".obj") continue;
        out.push_back({n.substr(0, n.size() - 4), dir + "/" + n, nullptr});
    }
    ::closedir(d);
    std::sort(out.begin(), out.end(),
              [](const Model &a, const Model &b) { return a.name < b.name; });
    return out;
}

// Load, returning null and leaving `why` set rather than throwing: a corpus
// sweep must survive one unreadable file.
inline std::shared_ptr<Mesh> load(const Model &m, std::string &why) {
    try {
        return std::make_shared<Mesh>(m.path);
    } catch (const std::exception &e) {
        why = e.what();
        return nullptr;
    }
}

// Keep only the models named on the command line, if any were.
inline std::vector<Model> select(std::vector<Model> all, const std::vector<std::string> &only) {
    if (only.empty()) return all;
    std::vector<Model> out;
    for (const Model &m : all)
        if (std::find(only.begin(), only.end(), m.name) != only.end()) out.push_back(m);
    return out;
}

} // namespace paper

#endif // __PAPER_CORPUS_HXX__
