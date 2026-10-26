#pragma once

#include <vector>
#include <cstdint>
#include <bit>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace cdbg {

class BitVector {
public:
    BitVector() = default;

    explicit BitVector(size_t num_bits) : num_bits_(num_bits) {
        size_t num_words = (num_bits + 63) / 64;
        words_.resize(num_words, 0);
    }

    size_t size() const noexcept { return num_bits_; }

    void set(size_t pos) {
        if (pos >= num_bits_) throw std::out_of_range("BitVector::set out of range");
        words_[pos / 64] |= (1ULL << (pos % 64));
    }

    void clear(size_t pos) {
        if (pos >= num_bits_) throw std::out_of_range("BitVector::clear out of range");
        words_[pos / 64] &= ~(1ULL << (pos % 64));
    }

    bool get(size_t pos) const {
        if (pos >= num_bits_) return false;
        return (words_[pos / 64] & (1ULL << (pos % 64))) != 0;
    }

    // Build 512-bit superblock ranks for O(1) rank1 queries
    void build_index() {
        size_t num_superblocks = (num_bits_ + 511) / 512;
        superblock_ranks_.resize(num_superblocks + 1, 0);

        uint64_t cumulative = 0;
        for (size_t sb = 0; sb < num_superblocks; ++sb) {
            superblock_ranks_[sb] = cumulative;
            size_t start_word = sb * 8;
            size_t end_word = std::min(start_word + 8, words_.size());
            for (size_t w = start_word; w < end_word; ++w) {
                cumulative += static_cast<uint64_t>(std::popcount(words_[w]));
            }
        }
        superblock_ranks_[num_superblocks] = cumulative;
        total_set_bits_ = cumulative;
    }

    // O(1) rank1: count number of 1-bits in range [0, pos)
    size_t rank1(size_t pos) const {
        if (pos == 0) return 0;
        if (pos > num_bits_) pos = num_bits_;

        size_t sb = pos / 512;
        uint64_t count = superblock_ranks_[sb];

        size_t start_word = sb * 8;
        size_t target_word = pos / 64;

        for (size_t w = start_word; w < target_word; ++w) {
            count += static_cast<uint64_t>(std::popcount(words_[w]));
        }

        size_t bit_offset = pos % 64;
        if (bit_offset > 0 && target_word < words_.size()) {
            uint64_t mask = (1ULL << bit_offset) - 1ULL;
            count += static_cast<uint64_t>(std::popcount(words_[target_word] & mask));
        }

        return count;
    }

    // select1(k): returns 0-based position of the k-th set bit (1 <= k <= total_set_bits)
    size_t select1(size_t k) const {
        if (k == 0 || k > total_set_bits_) {
            throw std::out_of_range("BitVector::select1 k out of range");
        }

        // Binary search on superblocks
        auto it = std::lower_bound(superblock_ranks_.begin(), superblock_ranks_.end(), k);
        size_t sb = (it == superblock_ranks_.begin()) ? 0 : static_cast<size_t>(std::distance(superblock_ranks_.begin(), it) - 1);

        uint64_t cur_rank = superblock_ranks_[sb];
        size_t start_word = sb * 8;
        size_t end_word = std::min(start_word + 8, words_.size());

        for (size_t w = start_word; w < end_word; ++w) {
            uint64_t pc = static_cast<uint64_t>(std::popcount(words_[w]));
            if (cur_rank + pc >= k) {
                // The bit is within words_[w]
                uint64_t word = words_[w];
                for (size_t b = 0; b < 64; ++b) {
                    if ((word & (1ULL << b)) != 0) {
                        cur_rank++;
                        if (cur_rank == k) {
                            return w * 64 + b;
                        }
                    }
                }
            }
            cur_rank += pc;
        }

        return num_bits_;
    }

    uint64_t total_set_bits() const noexcept { return total_set_bits_; }
    const std::vector<uint64_t>& words() const noexcept { return words_; }
    std::vector<uint64_t>& words() noexcept { return words_; }

    void save(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&num_bits_), sizeof(num_bits_));
        uint64_t nw = words_.size();
        out.write(reinterpret_cast<const char*>(&nw), sizeof(nw));
        if (nw > 0) {
            out.write(reinterpret_cast<const char*>(words_.data()), nw * sizeof(uint64_t));
        }
    }

    void load(std::istream& in) {
        in.read(reinterpret_cast<char*>(&num_bits_), sizeof(num_bits_));
        uint64_t nw = 0;
        in.read(reinterpret_cast<char*>(&nw), sizeof(nw));
        words_.resize(nw);
        if (nw > 0) {
            in.read(reinterpret_cast<char*>(words_.data()), nw * sizeof(uint64_t));
        }
        build_index();
    }

private:
    size_t num_bits_ = 0;
    std::vector<uint64_t> words_;
    std::vector<uint64_t> superblock_ranks_;
    uint64_t total_set_bits_ = 0;
};

} // namespace cdbg
