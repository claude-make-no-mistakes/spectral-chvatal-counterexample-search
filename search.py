"""Search rank-three and rank-four downsets, with numerical and exact output."""

from argparse import ArgumentParser
from pathlib import Path
from random import Random
import json

from certificates import make_certificate, verify
from theta import downset, layer, screen


def cyclic_orbits(n, rank):
    remaining = set(layer(n, rank))
    ceiling = (1 << n) - 1
    while remaining:
        first = min(remaining)
        orbit = {(first << k | first >> (n - k)) & ceiling for k in range(n)}
        remaining.difference_update(orbit)
        yield sorted(orbit)


def cases(seed, samples):
    random = Random(seed)
    for n in range(5, 11):
        for rank in (2, 3, 4):
            if n <= 8 or rank <= 3:
                yield "skeleton", n, rank, layer(n, rank)
    for _ in range(samples):
        n = random.choice((6, 7, 8, 9, 10))
        rank = random.choice((3, 3, 4))
        cyclic = random.choice((True, False))
        if cyclic:
            groups = list(cyclic_orbits(n, rank))
        else:
            groups = [[mask] for mask in layer(n, rank)]
        density = random.choice((0.15, 0.30, 0.55, 0.75, 0.90))
        selected = sorted(mask for group in groups if random.random() < density for mask in group)
        if not selected:
            continue
        if random.random() < 0.5:
            existing = set(downset(selected))
            selected.extend(mask for mask in layer(n, rank - 1) if mask not in existing)
        if len(downset(selected)) <= 225:
            yield "cyclic" if cyclic else "independent-facets", n, rank, selected


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=744103)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--certificate-dir", type=Path)
    args = parser.parse_args()
    if args.certificate_dir:
        args.certificate_dir.mkdir(parents=True, exist_ok=True)
    observed = set()
    print(json.dumps({"seed": args.seed, "requested_samples": args.samples}), flush=True)
    for family, n, rank, facets in cases(args.seed, args.samples):
        facets = tuple(sorted(facets))
        if facets in observed:
            continue
        observed.add(facets)
        result, residual = screen(facets)
        item = {key: value for key, value in result.items() if key not in ("primal", "slack")}
        item.update(family=family, n=n, rank=rank, facets=facets)
        item["gap"] = result["estimate"] - result["star"]
        item["gain_over_residual_star"] = item["gap"] + 1
        if args.certificate_dir and residual and str(result["status"]).startswith("solved"):
            direction = "lower" if item["gap"] > 1e-3 else "upper"
            if len(observed) <= 17 or item["gain_over_residual_star"] > 0.25:
                data = make_certificate(facets, residual, result, direction=direction)
                try:
                    item["certificate"] = verify(data)
                except ValueError:
                    item["certificate"] = None
                else:
                    path = args.certificate_dir / f"certificate-{len(observed):04d}.json"
                    path.write_text(json.dumps(data) + "\n", encoding="utf-8")
                    item["file"] = path.name
        print(json.dumps(item), flush=True)


if __name__ == "__main__":
    main()
