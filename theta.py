"""Screen hereditary set families for a strict Lovász theta gap."""

from itertools import combinations
from math import sqrt

import numpy as np
from scipy import sparse
import scs


def downset(facets):
    sets = set()
    for facet in facets:
        chunk = facet
        while chunk:
            sets.add(chunk)
            chunk = (chunk - 1) & facet
    return sorted(sets)


def star_size(sets):
    return max((sum(bool(mask & (1 << i)) for mask in sets)
                for i in range(max(sets, default=0).bit_length())), default=0)


def layer(n, k):
    return [sum(1 << i for i in row) for row in combinations(range(n), k)]


def theta(edges, *, tolerance=1e-6, iterations=7000):
    """Return a numerical screening estimate and primal/dual candidate matrices."""
    n = len(edges)
    rows, columns = np.triu_indices(n)
    diagonal = rows == columns
    packed_allowed = ~edges[rows, columns]
    positions = np.flatnonzero(packed_allowed)
    allowed_diagonal = diagonal[positions]
    scales = np.where(diagonal, 1.0, sqrt(2))
    count = len(positions)
    index = np.flatnonzero(allowed_diagonal)
    # Upper triangle in row order is the lower triangle in column order.
    A = sparse.csc_matrix((np.r_[np.ones(len(index)), -np.ones(count)],
                          (np.r_[np.zeros(len(index), dtype=int), positions + 1],
                           np.r_[index, np.arange(count)])), shape=(len(rows) + 1, count))
    b = np.zeros(len(rows) + 1)
    b[0] = 1
    c = -scales[positions]
    result = scs.solve({"A": A, "b": b, "c": c}, {"z": 1, "s": [n]},
                       verbose=False, eps_abs=tolerance, eps_rel=tolerance, max_iters=iterations)

    def unpack(triangle):
        matrix = np.zeros((n, n))
        matrix[rows, columns] = matrix[columns, rows] = triangle / scales
        return matrix

    point = np.zeros(len(rows))
    point[positions] = result["x"]
    info = result["info"]
    return {"estimate": float(-info["pobj"]), "dual": float(-info["dobj"]), "status": info["status"],
            "iterations": int(info["iter"]), "primal": unpack(point), "slack": unpack(result["y"][1:])}


def screen(facets, **options):
    sets = downset(facets)
    residual = [mask for mask in sets if mask.bit_count() >= 2]
    star = star_size(sets)
    if not residual:
        return {"star": star, "sets": len(sets), "residual": 0, "estimate": 0.0,
                "dual": 0.0, "status": "empty", "iterations": 0}, residual
    masks = np.array(residual, dtype=np.uint64)
    solution = theta((masks[:, None] & masks[None, :]) == 0, **options)
    return {"star": star, "sets": len(sets), "residual": len(residual), **solution}, residual
