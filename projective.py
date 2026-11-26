"""Replay the exceptional six-element family and its exact fractional gap."""

from argparse import ArgumentParser
from fractions import Fraction
from itertools import combinations, permutations
import json
from pathlib import Path

from certificates import verify
from small import closure, star_size


FACETS = (7, 11, 21, 26, 28, 38, 41, 44, 49, 50)
UNIVERSE = 63
CLIQUE_PARTITION = ((3, 12, 48), (5, 10, 16, 32), (7, 24), (11, 36), (21, 40, 2),
                    (26, 33, 4), (28, 34, 1), (38, 17, 8), (41, 20), (44, 18), (49, 6), (50, 9))


def permute(mask, permutation):
    return sum(1 << j for i, j in enumerate(permutation) if mask & (1 << i))


def orbit():
    return {sum(1 << a for a in closure(permute(f, p) for f in FACETS))
            for p in permutations(range(6))}


def partitions(remaining, allowed):
    if not remaining:
        yield ()
        return
    pivot = remaining & -remaining
    for a in allowed:
        if a & pivot and (a & remaining) == a:
            for rest in partitions(remaining ^ a, allowed):
                yield (a, *rest)


def verify_projective():
    sets = closure(FACETS)
    image_of_first_atom = {p[0] for p in permutations(range(6))
                           if {permute(a, p) for a in FACETS} == set(FACETS)}
    if image_of_first_atom != set(range(6)):
        raise ValueError("the base family must have a transitive automorphism group")
    if star_size(sets) != 11 or any((UNIVERSE ^ a) in sets for a in sets):
        raise ValueError("the listed facets did not give the exceptional downset")
    residual = sorted(a for a in sets if a.bit_count() >= 2)
    triangle = []
    for i, a in enumerate(residual):
        for j in range(i, len(residual)):
            b = residual[j]
            edge_weight = 0 if a & b else 2 if a.bit_count() == b.bit_count() == 2 else 4
            triangle.append(11 * (i == j) - 1 + edge_weight)
    theta_bound = verify({"facets": list(FACETS), "direction": "upper", "scale": 1, "triangle": triangle})

    # Fractions use denominator 3 in the packing and 9 in the covering.
    packing = {a: a.bit_count() - 1 for a in sets}
    all_partitions = list(partitions(UNIVERSE, sorted(sets)))
    if any(sum(packing[a] for a in p) > 3 for p in all_partitions):
        raise ValueError("the packing violates a clique constraint")
    cover = {a: 0 for a in sets}
    cost = 0
    for p in all_partitions:
        ranks = sorted(a.bit_count() for a in p)
        weight = 3 if ranks == [1, 2, 3] else 1 if ranks == [2, 2, 2] else 0
        cost += weight
        for a in p:
            cover[a] += weight
    if min(cover.values()) < 9 or Fraction(cost, 9) != Fraction(sum(packing.values()), 3):
        raise ValueError("the packing and covering certificates do not match")
    if sorted(a for part in CLIQUE_PARTITION for a in part) != sorted(sets) or any(
            a & b for part in CLIQUE_PARTITION for a, b in combinations(part, 2)):
        raise ValueError("the integral clique partition is invalid")
    return {**theta_bound, "fractional_independence": str(Fraction(cost, 9)),
            "minimum_clique_partition_size": len(CLIQUE_PARTITION),
            "partition_count": len(all_partitions), "labeled_forms": len(orbit())}


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--small6-log", type=Path)
    args = parser.parse_args()
    if args.small6_log:
        lines = args.small6_log.read_text(encoding="utf-8-sig").splitlines()
        cases = [line.split() for line in lines if line.startswith("undecided ") and len(line.split()) == 3]
        if (len(cases) != 12 or any(star != "11" for _, _, star in cases)
                or {int(word) for _, word, _ in cases} != orbit()
                or lines[-2:] != ["complement-free full-support downsets 1379272", "undecided 12"]):
            raise ValueError("the undecided cases are not exactly the projective-plane orbit")
    print(json.dumps(verify_projective()))
