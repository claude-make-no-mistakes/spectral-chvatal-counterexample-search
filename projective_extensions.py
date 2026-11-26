"""Screen seven-element extensions of the projective-plane downset."""

from argparse import ArgumentParser
from itertools import permutations
import json

from certificates import make_certificate, verify
from projective import FACETS, permute
from small import atom_union


def old_pairs():
    return [(1 << i) | (1 << j) for i in range(6) for j in range(i + 1, 6)]


def link_representatives():
    pairs = old_pairs()
    position = {a: i for i, a in enumerate(pairs)}
    actions = [tuple(1 << position[permute(a, p)] for a in pairs)
               for p in permutations(range(6)) if {permute(a, p) for a in FACETS} == set(FACETS)]
    seen = set()
    for state in range(1 << len(pairs)):
        if state in seen:
            continue
        indices = [i for i in range(len(pairs)) if state & (1 << i)]
        seen.update(sum(action[i] for i in indices) for action in actions)
        yield [pairs[i] for i in indices]


def extension_facets(link, attached):
    active = atom_union(link)
    return (*FACETS, *(64 | a for a in link),
            *(64 | (1 << i) for i in range(6) if (attached & ~active) & (1 << i)))


def extension_cases():
    for link in link_representatives():
        yield link, 63
    yield [a for a in old_pairs() if a < 32], 31


def main():
    from theta import screen

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--missing-join-only", action="store_true")
    args = parser.parse_args()
    for link, attached in extension_cases():
        if args.missing_join_only and attached == 63:
            continue
        facets = extension_facets(link, attached)
        result, residual = screen(facets)
        item = {k: v for k, v in result.items() if k not in ("primal", "slack")}
        item["link"] = link
        item["attached"] = attached
        item["gap"] = result["estimate"] - result["star"]
        if str(result["status"]).startswith("solved"):
            direction = "lower" if item["gap"] > 0.001 else "upper"
            data = make_certificate(facets, residual, result, direction=direction, scale=2000)
            try:
                item["exact"] = verify(data)
            except ValueError:
                pass
            else:
                item["certificate"] = data
        print(json.dumps(item), flush=True)


if __name__ == "__main__":
    main()
