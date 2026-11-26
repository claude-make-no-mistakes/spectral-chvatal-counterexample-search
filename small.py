"""Verify all small supports using complement pairs and residual clique covers."""

from argparse import ArgumentParser
from functools import cache
from itertools import combinations
import json


def antichains(n):
    universe = (1 << n) - 1
    incomparable = {a: sum(1 << b for b in range(1, universe + 1) if a & b not in (a, b))
                    for a in range(1, universe + 1)}

    def extend(remaining, facets):
        yield facets
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            a = bit.bit_length() - 1
            yield from extend(remaining & incomparable[a], (*facets, a))

    yield from extend((1 << (universe + 1)) - 2, ())


def closure(facets):
    sets = set()
    for facet in facets:
        subset = facet
        while subset:
            sets.add(subset)
            subset = (subset - 1) & facet
    return sets


def atom_union(sets):
    union = 0
    for mask in sets:
        union |= mask
    return union


def remove_complements(sets):
    """Return the remaining family and the number of removed pairs."""
    matched = 0
    while sets:
        universe = atom_union(sets)
        pairs = {a for a in sets if universe ^ a in sets}
        if not pairs:
            break
        matched += len(pairs) // 2
        sets = sets - pairs
    return sets, matched


def clique_cover_bound(sets):
    residual = sorted(a for a in sets if a.bit_count() >= 2)
    if atom_union(residual).bit_count() > 5:
        raise ValueError("the residual graph can have triangles on this support")
    neighbors = [sum(1 << j for j, b in enumerate(residual) if not (a & b)) for a in residual]

    @cache
    def matching(available):
        if not available:
            return 0
        first = available & -available
        other = available ^ first
        best = matching(other)
        choices = neighbors[first.bit_length() - 1] & other
        while choices:
            partner = choices & -choices
            choices ^= partner
            best = max(best, 1 + matching(other ^ partner))
        return best

    return len(residual) - matching((1 << len(residual)) - 1)


def star_size(sets):
    return max((sum(bool(a & (1 << i)) for a in sets)
                for i in range(max(sets, default=0).bit_length())), default=0)


def check(n):
    total = 0
    undecided = []
    for facets in antichains(n):
        total += 1
        family = closure(facets)
        reduced, pairs = remove_complements(family)
        star = star_size(reduced)
        if star + pairs != star_size(family):
            raise RuntimeError("removing complementary pairs must preserve the star deficit")
        bound = clique_cover_bound(reduced)
        if bound > star:
            undecided.append({"facets": facets, "reduced": sorted(reduced), "bound": bound, "star": star})
    return {"support_at_most": n, "nonempty_antichains": total - 1, "undecided": undecided}


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int, choices=range(1, 6), default=5, nargs="?")
    print(json.dumps(check(parser.parse_args().n)))
