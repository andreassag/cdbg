#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <algorithm>
#include <array>

namespace cdbg {

inline constexpr uint8_t base_to_2bit(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': case 'U': case 'u': return 3;
        default: return 4;
    }
}

inline constexpr char bit2_to_base(uint8_t b) noexcept {
    switch (b & 3) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

inline constexpr char complement_base(char c) noexcept {
    switch (c) {
        case 'A': return 'T';
        case 'a': return 't';
        case 'C': return 'G';
        case 'c': return 'g';
        case 'G': return 'C';
        case 'g': return 'c';
        case 'T': case 'U': return 'A';
        case 't': case 'u': return 'a';
        default: return 'N';
    }
}

inline std::string reverse_complement(std::string_view seq) {
    std::string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        rc.push_back(complement_base(*it));
    }
    return rc;
}

inline std::string canonical_kmer(std::string_view kmer) {
    std::string rc = reverse_complement(kmer);
    return (kmer < rc) ? std::string(kmer) : rc;
}

// 64-bit MurmurHash3 / SplitMix64 finalizer for string hashing
inline uint64_t hash64(std::string_view s, uint64_t seed = 0x9e3779b97f4a7c15ULL) noexcept {
    uint64_t h = seed ^ (static_cast<uint64_t>(s.size()) * 0x517cc1b727220a95ULL);
    for (char c : s) {
        h ^= static_cast<uint8_t>(c);
        h *= 0xbf58476d1ce4e5b9ULL;
        h ^= (h >> 31);
    }
    return h;
}

} // namespace cdbg
