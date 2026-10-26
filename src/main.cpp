#include "cdbg/fasta.hpp"
#include "cdbg/unitig_graph.hpp"
#include "cdbg/index_io.hpp"
#include <iostream>
#include <vector>
#include <string>

void print_usage() {
    std::cerr << "Usage: cdbg <command> [options]\n\n"
              << "Commands:\n"
              << "  build -k <int> -o <index.cdbg> <sample1.fa> [sample2.fa ...]\n"
              << "        Construct compacted colored de Bruijn graph and succinct filters\n"
              << "  query -i <index.cdbg> -q <query.fa>\n"
              << "        Query k-mers against the colored graph and succinct filters\n"
              << "  info -i <index.cdbg>\n"
              << "        Display graph topology, unitig count, and filter memory stats\n"
              << "  export-gfa -i <index.cdbg> -o <graph.gfa>\n"
              << "        Export unitigs and color annotations to GFA v1 format\n";
}

int cmd_build(int argc, char* argv[]) {
    uint32_t k = 31;
    std::string out_path;
    std::vector<std::string> sample_files;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-k" && i + 1 < argc) {
            k = static_cast<uint32_t>(std::stoul(argv[++i]));
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            out_path = argv[++i];
        } else if (arg[0] != '-') {
            sample_files.push_back(arg);
        }
    }

    if (sample_files.empty() || out_path.empty()) {
        std::cerr << "Error: At least one input FASTA and output path (-o) required.\n";
        print_usage();
        return 1;
    }

    std::vector<std::pair<std::string, std::string>> samples;
    for (const auto& file : sample_files) {
        auto recs = cdbg::read_fasta(file);
        for (const auto& r : recs) {
            samples.push_back({r.name, r.seq});
        }
    }

    std::cout << "[cdbg] Building graph for " << samples.size() << " sequences with k=" << k << "...\n";
    auto graph = cdbg::UnitigGraph::build(samples, k);

    std::cout << "[cdbg] Extracted " << graph.unitigs.size() << " unitigs, "
              << graph.colors.num_classes() << " color classes.\n";
    std::cout << "[cdbg] Bloom filter: " << graph.bloom_filter.size_in_bytes() << " bytes, "
              << "Xor filter: " << graph.xor_filter.size_in_bytes() << " bytes.\n";

    cdbg::save_cdbg(graph, out_path);
    std::cout << "[cdbg] Index saved to " << out_path << "\n";
    return 0;
}

int cmd_query(int argc, char* argv[]) {
    std::string index_path, query_path;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--index") && i + 1 < argc) {
            index_path = argv[++i];
        } else if ((arg == "-q" || arg == "--query") && i + 1 < argc) {
            query_path = argv[++i];
        }
    }

    if (index_path.empty() || query_path.empty()) {
        std::cerr << "Error: -i <index.cdbg> and -q <query.fa> required.\n";
        return 1;
    }

    auto graph = cdbg::load_cdbg(index_path);
    auto queries = cdbg::read_fasta(query_path);

    for (const auto& q : queries) {
        std::cout << "Query Sequence: " << q.name << " (" << q.seq.size() << " bp)\n";
        if (q.seq.size() < graph.k) {
            std::cout << "  (Sequence shorter than k=" << graph.k << ")\n";
            continue;
        }

        size_t total_kmers = q.seq.size() - graph.k + 1;
        size_t hits = 0;

        for (size_t i = 0; i < total_kmers; ++i) {
            std::string km = q.seq.substr(i, graph.k);
            auto res = graph.query_kmer(km);
            if (res.found) {
                hits++;
                if (hits <= 5) { // Print first 5 hits
                    std::cout << "  kmer[" << i << "] -> unitig " << res.unitig_id
                              << " at offset " << res.offset << " | colors: ";
                    for (size_t s = 0; s < res.present_samples.size(); ++s) {
                        if (s > 0) std::cout << ", ";
                        std::cout << res.present_samples[s];
                    }
                    std::cout << "\n";
                }
            }
        }
        std::cout << "  Found " << hits << " / " << total_kmers << " k-mers in graph ("
                  << (static_cast<double>(hits) * 100.0 / static_cast<double>(total_kmers)) << "%)\n";
    }

    return 0;
}

int cmd_info(int argc, char* argv[]) {
    std::string index_path;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--index") && i + 1 < argc) {
            index_path = argv[++i];
        }
    }
    if (index_path.empty()) {
        std::cerr << "Error: -i <index.cdbg> required.\n";
        return 1;
    }

    auto graph = cdbg::load_cdbg(index_path);
    size_t total_unitig_bp = 0;
    for (const auto& u : graph.unitigs) total_unitig_bp += u.seq.size();

    std::cout << "=== Colored Compacted De Bruijn Graph Info ===\n"
              << "k-mer size:            " << graph.k << "\n"
              << "Number of Unitigs:     " << graph.unitigs.size() << "\n"
              << "Total Unitig Bases:    " << total_unitig_bp << " bp\n"
              << "Number of Samples:     " << graph.colors.num_samples() << "\n"
              << "Color Equivalence Cls: " << graph.colors.num_classes() << "\n"
              << "Bloom Filter Size:     " << graph.bloom_filter.size_in_bytes() << " bytes\n"
              << "Xor Filter Size:       " << graph.xor_filter.size_in_bytes() << " bytes\n";
    return 0;
}

int cmd_export_gfa(int argc, char* argv[]) {
    std::string index_path, out_path;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--index") && i + 1 < argc) {
            index_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            out_path = argv[++i];
        }
    }
    if (index_path.empty() || out_path.empty()) {
        std::cerr << "Error: -i and -o required.\n";
        return 1;
    }

    auto graph = cdbg::load_cdbg(index_path);
    cdbg::export_gfa(graph, out_path);
    std::cout << "[cdbg] Exported GFA to " << out_path << "\n";
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }
    std::string cmd = argv[1];
    if (cmd == "build") return cmd_build(argc, argv);
    if (cmd == "query") return cmd_query(argc, argv);
    if (cmd == "info") return cmd_info(argc, argv);
    if (cmd == "export-gfa") return cmd_export_gfa(argc, argv);

    std::cerr << "Unknown command: " << cmd << "\n";
    print_usage();
    return 1;
}
