#pragma once

#include "bit_vector.hpp"
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

namespace cdbg {

class ColorMatrix {
public:
    ColorMatrix() = default;

    explicit ColorMatrix(std::vector<std::string> sample_names)
        : sample_names_(std::move(sample_names)) {}

    size_t num_samples() const noexcept { return sample_names_.size(); }
    const std::vector<std::string>& sample_names() const noexcept { return sample_names_; }

    // Registers a color bit-vector for a unitig, returns the equivalence class ID
    uint32_t get_or_create_class(const std::vector<bool>& sample_presence) {
        // Pack sample_presence into uint64_t or string key
        std::string key;
        for (bool b : sample_presence) key.push_back(b ? '1' : '0');

        auto it = class_map_.find(key);
        if (it != class_map_.end()) {
            return it->second;
        }

        uint32_t new_id = static_cast<uint32_t>(classes_.size());
        BitVector bv(sample_names_.size());
        for (size_t i = 0; i < sample_presence.size(); ++i) {
            if (sample_presence[i]) bv.set(i);
        }
        bv.build_index();
        classes_.push_back(std::move(bv));
        class_map_[key] = new_id;
        return new_id;
    }

    bool has_sample(uint32_t class_id, uint32_t sample_id) const {
        if (class_id >= classes_.size() || sample_id >= sample_names_.size()) {
            return false;
        }
        return classes_[class_id].get(sample_id);
    }

    size_t num_classes() const noexcept { return classes_.size(); }

    void save(std::ostream& out) const {
        uint32_t ns = static_cast<uint32_t>(sample_names_.size());
        out.write(reinterpret_cast<const char*>(&ns), sizeof(ns));
        for (const auto& name : sample_names_) {
            uint32_t len = static_cast<uint32_t>(name.size());
            out.write(reinterpret_cast<const char*>(&len), sizeof(len));
            out.write(name.data(), len);
        }

        uint32_t nc = static_cast<uint32_t>(classes_.size());
        out.write(reinterpret_cast<const char*>(&nc), sizeof(nc));
        for (const auto& bv : classes_) {
            bv.save(out);
        }
    }

    void load(std::istream& in) {
        uint32_t ns = 0;
        in.read(reinterpret_cast<char*>(&ns), sizeof(ns));
        sample_names_.resize(ns);
        for (size_t i = 0; i < ns; ++i) {
            uint32_t len = 0;
            in.read(reinterpret_cast<char*>(&len), sizeof(len));
            sample_names_[i].resize(len);
            in.read(&sample_names_[i][0], len);
        }

        uint32_t nc = 0;
        in.read(reinterpret_cast<char*>(&nc), sizeof(nc));
        classes_.resize(nc);
        class_map_.clear();
        for (size_t i = 0; i < nc; ++i) {
            classes_[i].load(in);
        }
    }

private:
    std::vector<std::string> sample_names_;
    std::vector<BitVector> classes_;
    std::unordered_map<std::string, uint32_t> class_map_;
};

} // namespace cdbg
