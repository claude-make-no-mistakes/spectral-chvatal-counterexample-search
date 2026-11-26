// Check complement-free downsets on a six-element support by exact clique covers.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {
using Word = std::uint64_t;
constexpr unsigned universe = 63;
std::array<Word, 64> subsets{}, compatible{}, disjoint{}, stars{};
std::array<std::vector<Word>, 64> triangles;
Word pairs = 0;
Word singletons = 0;
std::uint64_t seen = 0;
std::uint64_t undecided = 0;

constexpr Word bit(unsigned set) { return Word{1} << set; }
unsigned first(Word mask) { return static_cast<unsigned>(__builtin_ctzll(mask)); }
int size(Word mask) { return __builtin_popcountll(mask); }

bool cover(Word left, int saving) {
    if (saving <= 0) return true;
    const Word pair_vertices = left & pairs;
    if (saving > size(pair_vertices) || saving > 2 * size(left) / 3) return false;
    const unsigned a = first(pair_vertices);
    left ^= bit(a);
    for (const Word other_two : triangles[a]) {
        if ((left & other_two) == other_two && cover(left ^ other_two, saving - 2)) return true;
    }
    Word partners = disjoint[a] & left;
    while (partners) {
        const Word partner = bit(first(partners));
        partners ^= partner;
        if (cover(left ^ partner, saving - 1)) return true;
    }
    return cover(left, saving);
}

void verify(Word family, unsigned support) {
    if (support != universe) return;
    ++seen;
    int star = 0;
    for (unsigned i = 0; i < 6; ++i) {
        const int count = size(family & stars[i]);
        if (count > star) star = count;
    }
    const Word residual = family & ~singletons;
    if (!cover(residual, size(residual) - star)) {
        ++undecided;
        std::cout << "undecided " << family << ' ' << star << '\n';
    }
}

void enumerate(Word available, Word family, unsigned support) {
    verify(family, support);
    while (available) {
        const unsigned a = first(available);
        available ^= bit(a);
        enumerate(available & compatible[a], family | subsets[a], support | a);
    }
}
}

int main() {
    for (unsigned a = 1; a < universe; ++a) {
        if (size(a) == 1) singletons |= bit(a);
        if (size(a) == 2) pairs |= bit(a);
        for (unsigned i = 0; i < 6; ++i) {
            if (a & (1u << i)) stars[i] |= bit(a);
        }
        for (unsigned b = 1; b < universe; ++b) {
            const unsigned intersection = a & b;
            if (intersection == b) subsets[a] |= bit(b);
            if (intersection != a && intersection != b && (a | b) != universe) compatible[a] |= bit(b);
            if (!intersection) disjoint[a] |= bit(b);
        }
    }
    for (unsigned a = 1; a < universe; ++a) {
        if (size(a) != 2) continue;
        for (unsigned b = 1; b < universe; ++b) {
            if (size(b) != 2 || (a & b)) continue;
            const unsigned c = universe ^ (a | b);
            if (b < c) triangles[a].push_back(bit(b) | bit(c));
        }
    }
    enumerate(~(bit(0) | bit(universe)), 0, 0);
    std::cout << "complement-free full-support downsets " << seen << "\nundecided " << undecided << '\n';
}
