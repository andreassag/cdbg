#pragma once

#include <vector>
#include <cstdint>
#include <array>
#include <bit>
#include <cmath>
#include <stdexcept>
#include <cstring>
#include <iostream>

namespace cdbg {

// Fast 64-bit integer mix
inline uint64_t splitmix64(uint64_t x) noexcept {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

// -------------------------------------------------------------
// Cache-Line Blocked Bloom Filter (512 bits = 64 bytes per block)
// -------------------------------------------------------------
class BlockedBloomFilter {
public:
    struct alignas(64) Block {
        uint64_t words[8] = {0};
    };

    BlockedBloomFilter() = default;

    explicit BlockedBloomFilter(size_t expected_elements, double target_fpr = 0.01) {
        // Optimal bits: m = -n * ln(p) / (ln2)^2 ~= 10 bits/elem for 1% FPR
        double bits_per_elem = -std::log(target_fpr) / (0.693147 * 0.693147);
        if (bits_per_elem < 8.0) bits_per_elem = 8.0;
        size_t total_bits = static_cast<size_t>(static_cast<double>(expected_elements) * bits_per_elem);
        size_t num_blocks = (total_bits + 511) / 512;
        if (num_blocks == 0) num_blocks = 1;
        blocks_.resize(num_blocks);
    }

    void insert(uint64_t hash) {
        if (blocks_.empty()) return;
        size_t block_idx = (hash >> 32) % blocks_.size();
        Block& blk = blocks_[block_idx];

        uint64_t h = hash;
        for (int i = 0; i < 4; ++i) {
            h = splitmix64(h);
            size_t bit_idx = static_cast<size_t>(h % 512);
            blk.words[bit_idx / 64] |= (1ULL << (bit_idx % 64));
        }
    }

    bool contains(uint64_t hash) const {
        if (blocks_.empty()) return false;
        size_t block_idx = (hash >> 32) % blocks_.size();
        const Block& blk = blocks_[block_idx];

        uint64_t h = hash;
        for (int i = 0; i < 4; ++i) {
            h = splitmix64(h);
            size_t bit_idx = static_cast<size_t>(h % 512);
            if ((blk.words[bit_idx / 64] & (1ULL << (bit_idx % 64))) == 0) {
                return false;
            }
        }
        return true;
    }

    size_t num_blocks() const noexcept { return blocks_.size(); }
    size_t size_in_bytes() const noexcept { return blocks_.size() * sizeof(Block); }

    void save(std::ostream& out) const {
        uint64_t nb = blocks_.size();
        out.write(reinterpret_cast<const char*>(&nb), sizeof(nb));
        if (nb > 0) {
            out.write(reinterpret_cast<const char*>(blocks_.data()), nb * sizeof(Block));
        }
    }

    void load(std::istream& in) {
        uint64_t nb = 0;
        in.read(reinterpret_cast<char*>(&nb), sizeof(nb));
        blocks_.resize(nb);
        if (nb > 0) {
            in.read(reinterpret_cast<char*>(blocks_.data()), nb * sizeof(Block));
        }
    }

private:
    std::vector<Block> blocks_;
};

// -------------------------------------------------------------
// Xor Filter 8-bit (Graf & Lemire 2020)
// Zero false negatives, ~0.39% FPR, ~9.84 bits per key
// -------------------------------------------------------------
class XorFilter8 {
public:
    XorFilter8() = default;

    explicit XorFilter8(const std::vector<uint64_t>& keys) {
        build(keys);
    }

    bool build(const std::vector<uint64_t>& keys) {
        if (keys.empty()) {
            fingerprints_.clear();
            size_ = 0;
            return true;
        }

        size_t n = keys.size();
        // Graf & Lemire factor: 1.23 * n + 32
        size_t capacity = static_cast<size_t>(1.23 * static_cast<double>(n)) + 32;
        capacity_ = (capacity / 3) * 3; // Ensure divisible by 3
        segment_len_ = capacity_ / 3;

        fingerprints_.assign(capacity_, 0);

        uint64_t seed = 0x85432901ULL;
        for (int trial = 0; trial < 100; ++trial) {
            seed = splitmix64(seed + static_cast<uint64_t>(trial));
            if (try_peel(keys, seed)) {
                seed_ = seed;
                size_ = n;
                return true;
            }
        }

        throw std::runtime_error("XorFilter8: Failed to peel hypergraph after 100 trials");
    }

    bool contains(uint64_t key) const {
        if (capacity_ == 0) return false;

        auto [h0, h1, h2] = get_hashes(key, seed_);
        uint8_t fp = get_fingerprint(key);
        return (fingerprints_[h0] ^ fingerprints_[h1] ^ fingerprints_[h2]) == fp;
    }

    size_t size() const noexcept { return size_; }
    size_t capacity() const noexcept { return capacity_; }
    size_t size_in_bytes() const noexcept { return fingerprints_.size(); }

    void save(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&size_), sizeof(size_));
        out.write(reinterpret_cast<const char*>(&capacity_), sizeof(capacity_));
        out.write(reinterpret_cast<const char*>(&segment_len_), sizeof(segment_len_));
        out.write(reinterpret_cast<const char*>(&seed_), sizeof(seed_));
        uint64_t fp_size = fingerprints_.size();
        out.write(reinterpret_cast<const char*>(&fp_size), sizeof(fp_size));
        if (fp_size > 0) {
            out.write(reinterpret_cast<const char*>(fingerprints_.data()), fp_size);
        }
    }

    void load(std::istream& in) {
        in.read(reinterpret_cast<char*>(&size_), sizeof(size_));
        in.read(reinterpret_cast<char*>(&capacity_), sizeof(capacity_));
        in.read(reinterpret_cast<char*>(&segment_len_), sizeof(segment_len_));
        in.read(reinterpret_cast<char*>(&seed_), sizeof(seed_));
        uint64_t fp_size = 0;
        in.read(reinterpret_cast<char*>(&fp_size), sizeof(fp_size));
        fingerprints_.resize(fp_size);
        if (fp_size > 0) {
            in.read(reinterpret_cast<char*>(fingerprints_.data()), fp_size);
        }
    }

private:
    size_t size_ = 0;
    size_t capacity_ = 0;
    size_t segment_len_ = 0;
    uint64_t seed_ = 0;
    std::vector<uint8_t> fingerprints_;

    std::array<size_t, 3> get_hashes(uint64_t key, uint64_t seed) const noexcept {
        uint64_t h0 = splitmix64(key ^ seed);
        uint64_t h1 = splitmix64(h0);
        uint64_t h2 = splitmix64(h1);

        size_t p0 = static_cast<size_t>(h0 % segment_len_);
        size_t p1 = static_cast<size_t>(segment_len_ + (h1 % segment_len_));
        size_t p2 = static_cast<size_t>(2 * segment_len_ + (h2 % segment_len_));
        return {p0, p1, p2};
    }

    uint8_t get_fingerprint(uint64_t key) const noexcept {
        uint64_t h = splitmix64(key ^ 0xa0761d6478bd642fULL);
        uint8_t fp = static_cast<uint8_t>((h ^ (h >> 32)) & 0xFF);
        return (fp == 0) ? 1 : fp; // ensure non-zero
    }

    struct KeyEdge {
        uint64_t key;
        size_t h[3];
    };

    bool try_peel(const std::vector<uint64_t>& keys, uint64_t seed) {
        size_t n = keys.size();
        std::vector<uint32_t> degree(capacity_, 0);
        std::vector<uint64_t> xor_keys(capacity_, 0);

        for (uint64_t k : keys) {
            auto [h0, h1, h2] = get_hashes(k, seed);
            degree[h0]++; degree[h1]++; degree[h2]++;
            xor_keys[h0] ^= k;
            xor_keys[h1] ^= k;
            xor_keys[h2] ^= k;
        }

        std::vector<size_t> queue;
        queue.reserve(capacity_);
        for (size_t i = 0; i < capacity_; ++i) {
            if (degree[i] == 1) {
                queue.push_back(i);
            }
        }

        struct StackEntry {
            size_t vertex;
            uint64_t key;
        };
        std::vector<StackEntry> stack;
        stack.reserve(n);

        size_t head = 0;
        while (head < queue.size()) {
            size_t v = queue[head++];
            if (degree[v] == 0) continue;

            uint64_t k = xor_keys[v];
            stack.push_back({v, k});

            auto [h0, h1, h2] = get_hashes(k, seed);
            for (size_t h : {h0, h1, h2}) {
                degree[h]--;
                xor_keys[h] ^= k;
                if (degree[h] == 1) {
                    queue.push_back(h);
                }
            }
        }

        if (stack.size() < n) {
            // Cycle prevented full peeling
            return false;
        }

        // Backtrack stack and solve for fingerprints
        fingerprints_.assign(capacity_, 0);
        while (!stack.empty()) {
            auto [v, k] = stack.back();
            stack.pop_back();

            auto [h0, h1, h2] = get_hashes(k, seed);
            uint8_t fp = get_fingerprint(k);

            // fingerprints[h0] ^ fingerprints[h1] ^ fingerprints[h2] = fp
            // v is one of {h0, h1, h2}
            uint8_t current_xor = 0;
            for (size_t h : {h0, h1, h2}) {
                if (h != v) {
                    current_xor ^= fingerprints_[h];
                }
            }
            fingerprints_[v] = fp ^ current_xor;
        }

        return true;
    }
};

} // namespace cdbg
