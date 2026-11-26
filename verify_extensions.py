"""Replay certificates for all rank-three one-point extensions of the projective plane."""

from argparse import ArgumentParser
import gzip
import json
from pathlib import Path

from certificates import verify
from projective_extensions import extension_cases, extension_facets


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    count = 0
    with gzip.open(args.archive, "rt", encoding="utf-8") as stream:
        for (link, attached), row in zip(extension_cases(), map(json.loads, stream), strict=True):
            certificate = row["certificate"]
            if row["link"] != link or row.get("attached", 63) != attached:
                raise ValueError("certificate rows did not match the link enumeration")
            if certificate["facets"] != list(extension_facets(link, attached)) or certificate["direction"] != "upper":
                raise ValueError("wrong family or certificate direction in the archive")
            verify(certificate)
            count += 1
    print(json.dumps({"upper_certificates_verified": count}))
