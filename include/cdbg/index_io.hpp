#pragma once

#include "unitig_graph.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <iomanip>

namespace cdbg {

inline void save_cdbg(const UnitigGraph& graph, const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    if (!out.is_open()) {
        throw std::runtime_error("Could not write .cdbg file: " + path);
    }

    const char magic[4] = {'C', 'D', 'B', 'G'};
    uint16_t version = 1;
    out.write(magic, 4);
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));

    out.write(reinterpret_cast<const char*>(&graph.k), sizeof(graph.k));

    // Unitigs
    uint64_t num_unitigs = graph.unitigs.size();
    out.write(reinterpret_cast<const char*>(&num_unitigs), sizeof(num_unitigs));
    for (const auto& u : graph.unitigs) {
        out.write(reinterpret_cast<const char*>(&u.id), sizeof(u.id));
        out.write(reinterpret_cast<const char*>(&u.color_class), sizeof(u.color_class));
        uint32_t len = static_cast<uint32_t>(u.seq.size());
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        out.write(u.seq.data(), len);

        uint32_t num_edges = static_cast<uint32_t>(u.out_edges.size());
        out.write(reinterpret_cast<const char*>(&num_edges), sizeof(num_edges));
        if (num_edges > 0) {
            out.write(reinterpret_cast<const char*>(u.out_edges.data()), num_edges * sizeof(uint32_t));
        }
    }

    // Color matrix
    graph.colors.save(out);

    // Filters
    graph.xor_filter.save(out);
    graph.bloom_filter.save(out);
}

inline UnitigGraph load_cdbg(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("Could not read .cdbg file: " + path);
    }

    char magic[4];
    in.read(magic, 4);
    if (std::memcmp(magic, "CDBG", 4) != 0) {
        throw std::runtime_error("Invalid .cdbg magic header: " + path);
    }

    uint16_t version = 0;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != 1) {
        throw std::runtime_error("Unsupported .cdbg version: " + std::to_string(version));
    }

    UnitigGraph graph;
    in.read(reinterpret_cast<char*>(&graph.k), sizeof(graph.k));

    uint64_t num_unitigs = 0;
    in.read(reinterpret_cast<char*>(&num_unitigs), sizeof(num_unitigs));
    graph.unitigs.resize(num_unitigs);

    for (size_t i = 0; i < num_unitigs; ++i) {
        auto& u = graph.unitigs[i];
        in.read(reinterpret_cast<char*>(&u.id), sizeof(u.id));
        in.read(reinterpret_cast<char*>(&u.color_class), sizeof(u.color_class));
        uint32_t len = 0;
        in.read(reinterpret_cast<char*>(&len), sizeof(len));
        u.seq.resize(len);
        in.read(&u.seq[0], len);

        uint32_t num_edges = 0;
        in.read(reinterpret_cast<char*>(&num_edges), sizeof(num_edges));
        u.out_edges.resize(num_edges);
        if (num_edges > 0) {
            in.read(reinterpret_cast<char*>(u.out_edges.data()), num_edges * sizeof(uint32_t));
        }
    }

    graph.colors.load(in);
    graph.xor_filter.load(in);
    graph.bloom_filter.load(in);

    return graph;
}

inline void export_gfa(const UnitigGraph& graph, const std::string& path) {
    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Could not write GFA file: " + path);
    }

    out << "H\tVN:Z:1.0\n";

    for (const auto& u : graph.unitigs) {
        out << "S\tutg" << std::setw(6) << std::setfill('0') << u.id << "\t"
            << u.seq << "\t"
            << "LN:i:" << u.seq.size() << "\t"
            << "colors:Z:";

        bool first = true;
        for (size_t s = 0; s < graph.colors.num_samples(); ++s) {
            if (graph.colors.has_sample(u.color_class, static_cast<uint32_t>(s))) {
                if (!first) out << ",";
                out << graph.colors.sample_names()[s];
                first = false;
            }
        }
        out << "\n";
    }

    // Export links
    for (const auto& u : graph.unitigs) {
        for (uint32_t tgt : u.out_edges) {
            out << "L\tutg" << std::setw(6) << std::setfill('0') << u.id << "\t+\t"
                << "utg" << std::setw(6) << std::setfill('0') << tgt << "\t+\t"
                << (graph.k - 1) << "M\n";
        }
    }
}

} // namespace cdbg
