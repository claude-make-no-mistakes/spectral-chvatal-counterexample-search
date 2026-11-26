# Spectral Chvátal on six elements

Every downset on at most six elements satisfies the Lovász theta equality in Conjecture H. Exactly one family, up to relabeling, violates the stronger fractional independence bound in this range. The face family of the six-vertex projective plane has largest star and theta value both equal to 11, fractional independence $35/3$, and minimum clique partition size 12. Six is the least support size for a fractional gap. [Paper](paper.pdf) · [Proof source](paper.tex) · [Conjecture H and the earlier seven-element fractional example](https://arxiv.org/html/2609.28404v3#S5)

The unrestricted Conjecture H remains open. The search reduces it to families with no complementary pairs and an SDP on sets of size at least two. A strict primal gap can be checked by integer arithmetic after rounding and a small diagonal shift. An exact rational dual upper bound also uses only integer arithmetic in the checker.

The same equality holds for every rank-three extension of the exceptional family by a seventh element. The archive contains exact theta duals for the 688 link graph orbits with all six new pairs present and the one extra missing-pair case where the monotonicity argument needs a separate certificate.

## Exact replay

The finite proof uses Python 3.10 or newer and a GCC-compatible C++17 compiler. It requires no numerical solver. The five-element check covers 7,579 labeled nontrivial downsets. The six-element check covers all 1,379,272 labeled full-support complement-free families. The twelve entries undecided by a clique partition constitute one orbit. The verifier checks that orbit, then establishes theta and both clique partition values on its representative.

```sh
python -B small.py 5
g++ -O3 -std=c++17 -Wall -Wextra -Wpedantic small6.cpp -o small6
./small6 > small6.log
python -B projective.py --small6-log small6.log
python -B verify_extensions.py projective-extensions.jsonl.gz
```

## Search larger examples

Numerical screening uses the pinned dependencies and keeps the SDP estimates separate from exact matrix certificates.

```sh
python -m pip install -r requirements-screening.txt
python -B search.py --seed 744103 --samples 130 --certificate-dir certificates
python -B projective_extensions.py > projective-extensions.jsonl
```
