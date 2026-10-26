#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

namespace cdbg {

struct FastaRecord {
    std::string name;
    std::string comment;
    std::string seq;
};

inline std::vector<FastaRecord> read_fasta(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open FASTA file: " + path);
    }

    std::vector<FastaRecord> records;
    std::string line;
    FastaRecord cur;
    bool in_record = false;

    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (in_record) {
                records.push_back(std::move(cur));
                cur = FastaRecord{};
            }
            in_record = true;
            size_t space_pos = line.find_first_of(" \t");
            if (space_pos != std::string::npos) {
                cur.name = line.substr(1, space_pos - 1);
                cur.comment = line.substr(space_pos + 1);
            } else {
                cur.name = line.substr(1);
            }
        } else if (in_record) {
            for (char c : line) {
                if (c > ' ') {
                    cur.seq.push_back(c);
                }
            }
        }
    }

    if (in_record) {
        records.push_back(std::move(cur));
    }
    return records;
}

} // namespace cdbg
