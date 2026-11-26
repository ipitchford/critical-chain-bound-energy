# SPDX-License-Identifier: MIT
"""Compare regenerated results with the shipped results.

Usage: python3 compare_results.py SHIPPED_DIR NEW_DIR
Policy (also stated in README.md, ENVIRONMENT.txt and REPLAY_RECEIPT.md):
- five_site_certificate.txt and ising_boundary.csv must be byte-identical (exact interval and integer outputs);
- every other file must have the same rows, columns and non-numeric fields; numeric fields must agree to a
  relative tolerance of 1e-6, except values below 1e-6 in magnitude (error diagnostics and small differences), which
  must agree to an absolute tolerance of 1e-7;
- finite_size.csv is compared at relative 1e-3: its L = 4096 one-body column for the cubic chain sums level
  energies of order 1e-10 from a dense eigensolver and changed by about 2e-4 (relative) between Linux runners.
Exits 1 on any difference. Uses no assert, so it also fails under python -O.
"""
import csv, sys
from pathlib import Path

EXACT = {'five_site_certificate.txt', 'ising_boundary.csv'}
REL, SMALL, ABS = 1e-6, 1e-6, 1e-7
LOOSE = {'finite_size.csv': 1e-3}


def rows(p):
    with open(p, newline='') as f:
        return list(csv.reader(f))


def close(u, v, rel=REL):
    if u == v:
        return True
    try:
        a, b = float(u), float(v)
    except ValueError:
        return False
    m = max(abs(a), abs(b))
    if m < SMALL:
        return abs(a - b) <= ABS
    return abs(a - b) <= rel * m


def main():
    old, new = Path(sys.argv[1]), Path(sys.argv[2])
    names = sorted(p.name for p in old.iterdir() if p.is_file())
    if names != sorted(p.name for p in new.iterdir() if p.is_file()):
        print('FAIL file sets differ'); sys.exit(1)
    bad = []
    for n in names:
        if n in EXACT:
            if (old / n).read_bytes() != (new / n).read_bytes():
                bad.append(n + ' not byte-identical')
            continue
        a, b = rows(old / n), rows(new / n)
        if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
            bad.append(n + ' shape differs'); continue
        for i, (x, y) in enumerate(zip(a, b)):
            for u, v in zip(x, y):
                if not close(u, v, LOOSE.get(n, REL)):
                    bad.append('%s row %d: %s vs %s' % (n, i, u, v))
    if bad:
        print('FAIL', bad[:20]); sys.exit(1)
    print('PASS %d files: exact files byte-identical; others within rel 1e-6 (abs 1e-7 below 1e-6; finite_size.csv rel 1e-3)' % len(names))


main()
