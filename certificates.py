"""Export and verify integer matrix certificates with no numerical dependencies."""

from argparse import ArgumentParser
from fractions import Fraction
from pathlib import Path
import json


def positive_definite(matrix):
    """Use fraction-free elimination to test strictly positive leading minors."""
    work = [row.copy() for row in matrix]
    previous = 1
    for k, row in enumerate(work):
        pivot = row[k]
        if pivot <= 0:
            return False
        for i in range(k + 1, len(work)):
            for j in range(i, len(work)):
                numerator = work[i][j] * pivot - row[i] * row[j]
                value, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("non-exact Bareiss division")
                work[i][j] = value
        previous = pivot
    return True


def make_certificate(facets, residual, result, *, direction, scale=100_000):
    source = result["slack" if direction == "upper" else "primal"]
    matrix = []
    for i, a in enumerate(residual):
        row = []
        for j, b in enumerate(residual):
            if direction == "upper" and a & b:
                value = ((result["star"] if i == j else 0) - 1) * scale
            elif direction == "lower" and not (a & b):
                value = 0
            else:
                value = int(round(scale * source[i, j]))
            row.append(value)
        matrix.append(row)
    if direction == "lower" and matrix:
        # A diagonal shift preserves the disjoint zeros. Try it only while
        # the normalization by trace can still leave a strict objective gap.
        n = len(matrix)
        trace = sum(matrix[i][i] for i in range(n))
        surplus = sum(map(sum, matrix)) - result["star"] * trace
        shift = 1
        while surplus > (result["star"] - 1) * n * shift and not positive_definite(matrix):
            for i in range(n):
                matrix[i][i] += shift
            surplus -= (result["star"] - 1) * n * shift
            shift *= 2
    data = {"facets": list(facets), "direction": direction, "scale": scale,
            "triangle": [matrix[i][j] for i in range(len(matrix)) for j in range(i, len(matrix))]}
    return data


def verify(data):
    facets = data["facets"]
    if (not facets or any(type(facet) is not int or facet < 1 for facet in facets)
            or len(set(facets)) != len(facets)):
        raise ValueError("facets must be distinct positive integer bitmasks")
    # Build the downset by taking one-step subsets until closure.
    stack = facets.copy()
    sets = set()
    while stack:
        mask = stack.pop()
        if mask in sets:
            continue
        sets.add(mask)
        subset_bits = mask
        while subset_bits:
            bit = subset_bits & -subset_bits
            subset = mask ^ bit
            if subset:
                stack.append(subset)
            subset_bits -= bit
    stars = [sum(bool((mask >> element) & 1) for mask in sets)
             for element in range(max(sets).bit_length())]
    star = max(stars)
    residual = sorted(mask for mask in sets if mask.bit_count() >= 2)
    size = len(residual)
    scale = data["scale"]
    upper = data["direction"] == "upper"
    if data["direction"] not in ("upper", "lower") or type(scale) is not int or scale < 1:
        raise ValueError("invalid certificate type or scale")
    triangle = data["triangle"]
    if len(triangle) != size * (size + 1) // 2 or any(type(value) is not int for value in triangle):
        raise ValueError("triangle must contain exactly the upper-triangular integer entries")
    values = iter(triangle)
    matrix = [[0] * size for _ in residual]
    for i, a in enumerate(residual):
        for j in range(i, size):
            b = residual[j]
            value = next(values)
            if upper and a & b and value != ((star if i == j else 0) - 1) * scale:
                raise ValueError("incompatible dual entry")
            if not upper and not (a & b) and value != 0:
                raise ValueError("incompatible primal entry")
            matrix[i][j] = matrix[j][i] = value
    if not positive_definite(matrix):
        raise ValueError("certificate matrix is not positive definite")
    if upper:
        return {"star": star, "residual": size, "upper_bound": star}
    trace = sum(matrix[i][i] for i in range(size))
    if not trace:
        raise ValueError("the primal certificate is empty")
    total = sum(map(sum, matrix))
    if total <= star * trace:
        raise ValueError("the primal objective does not exceed the largest star")
    return {"star": star, "residual": size, "lower_bound": str(Fraction(total, trace))}


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path, nargs="+")
    args = parser.parse_args()
    for path in args.certificate:
        print(json.dumps({"file": str(path), **verify(json.loads(path.read_text(encoding="utf-8")))}))
