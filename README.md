# The energy an isolated half of a critical chain cannot give up

*Conformal amplitudes, dispersion, and the cost of conserving particle number*

Version 0.1.0-candidate, 26 November 2026. DOI [10.5281/zenodo.24553557](https://doi.org/10.5281/zenodo.24553557)
(all versions: [10.5281/zenodo.24553556](https://doi.org/10.5281/zenodo.24553556)). Author: Anonymous.
Release page: https://evidencepress.org/releases/critical-chain-bound-energy/

**Status: unrefereed research candidate.** The results are written proofs with computational checks. They were
checked by internal AI-assisted review, proof audits and replays only; see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
This is not peer review, independent reproduction or formal verification.

## What this is about

Take a long chain of quantum particles in its lowest-energy state and cut it in half. Each half is entangled with the
other, so on its own it is not in its own lowest-energy state: it holds some extra energy. Operations on that half
alone can pull most of the extra energy out as useful work (the *ergotropy*), but not all of it. The part that can
never be pulled out is the *bound energy* Q_L of a half of length L. It measures how much energy entanglement locks
away.

Mula and colleagues (Phys. Rev. B 107, 075116, 2023) computed the bound energy for the simplest critical chain, free
fermions hopping between neighbouring sites, and found that it scales like the square of the entanglement entropy S
divided by the length, N Q ≈ (6/π) S² for a chain of N = 2L sites. They conjectured the same law for critical chains
in general. Later work (Mitra and Srivastava, New J. Phys. 27, 084601, 2025) saw the law beyond conformal symmetry in
a long-range Kitaev chain and a deviation for all-to-all pairing, but no form of the conjecture had been proved.

## What the paper shows

- **The conformal law, with its exact coefficient.** For a regulated boundary conformal field theory,
  Q_L ~ (3v / (2πc)) S_L² / L, where v is the speed of excitations and c the central charge. The same law is proved
  directly for the open hopping chain and the critical Ising chain: slopes 3/π and 6/π in units of L (6/π and 12/π in
  units of N). At equal speed the Ising slope is twice the hopping slope, so the conjectured number is not universal.
- **One ground state, many bound energies.** A whole family of free-fermion chains shares exactly the hopping chain's
  critical ground state, yet their bound energies scale with different powers of L: Q_L ~ D S_L^(m+1) / L^m, where m
  is the order of the flattest zero of the cut-off half's own one-particle energy band. The familiar S²/L law holds
  exactly when the band crosses the Fermi level with nonzero speed and touches zero nowhere else. So the broadest form
  of the conjecture ("all critical states") is false. (Hartmann, Kavanagh and Vandoren had shown that odd Lifshitz
  dispersions leave the ground-state entanglement unchanged; the new part is the exact bound-energy law.)
- **Stability.** Adding a linear term to a flat dispersion restores the S²/L law through an explicit crossover
  function.
- **What you are allowed to do matters.** In one range-three chain, unrestricted operations, particle-number-conserving
  operations and number-conserving one-body rotations leave residual energies of different orders:
  (log L)³/L², (log L)/L and (log L)²/L. The number-conserving remainder is already at least about 20 times the
  unrestricted one at L = 16 and about 1,360 times at L = 4096 (Table 4).

**Idea of the proofs.** A rank-matching lemma shows that, to leading order, the bound energy equals the energy at which
the half's own thermal state has entropy equal to its entanglement entropy. The work is in counting how many energy
levels lie below a given energy and in controlling the rare, improbable Schmidt ranks.

**Limits.** The lattice proofs cover free-fermion chains; interacting chains such as the Heisenberg chain remain open.
They use one external input, the known boundary Rényi asymptotics of the open hopping chain. The optimal operations are
idealised (no limit on circuit depth, locality or duration). All coefficients are limits as L grows, and some are
reached only logarithmically slowly: for the cubic chain the one-body penalty is about 1.096 at L = 4096 against a
limit of about 2.148.

## Files

| File | Contents |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper and its LaTeX source (figures in [figures/](figures/)) |
| [code/check_passive_fermi.py](code/check_passive_fermi.py) | Complete Fock-space ordering and direct many-body checks for odd-power chains → results/passive_fermi.csv |
| [code/check_asymmetric_touch.py](code/check_asymmetric_touch.py) | Range-three chain: unrestricted, charge-conserving and one-body residuals → results/asymmetric_touch.csv |
| [code/certify_five_sites.py](code/certify_five_sites.py) | Interval-arithmetic certificate (mpmath) of the ordering inversion at five sites → results/five_site_certificate.txt |
| [code/check_boundary_characters.py](code/check_boundary_characters.py) | Ising annulus spectra with exact integer multiplicities → results/ising_boundary.csv |
| [code/reproduce_ising.py](code/reproduce_ising.py), [code/ising_mesoscopic.py](code/ising_mesoscopic.py) | Ising spin and Majorana calculations → results/ising_checks.csv, results/ising_mesoscopic.csv |
| [code/crossover.py](code/crossover.py), [code/dispersion.py](code/dispersion.py) | Crossover quadrature and Figures 1–2 → results/crossover.csv, figures/ |
| [code/general_dispersion_check.py](code/general_dispersion_check.py) | Pressure constants of the general-dispersion theorem and the constants quoted in the text → results/general_dispersion.csv |
| [code/finite_size.py](code/finite_size.py) | Finite-size table up to L = 4096 (Table 4) → results/finite_size.csv |
| [tools/compare_results.py](tools/compare_results.py) | Compares regenerated results with the shipped ones under the policy below |
| [refs.bib](refs.bib) | Bibliography (16 entries) |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits |
| [RESEARCH_GATES.json](RESEARCH_GATES.json) | Prior-work map, sources read, limitation triage |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | Public summary of the internal review, without scores |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Release replay and environment |
| [CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256) | Citation, licences, SHA-256 sums |

## Reproduce

From the package root, in a copy (the scripts rewrite `results/` and `figures/` next to `code/`), about a minute on
one CPU:

```sh
export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
U="uv run --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with mpmath==1.3.0 --with matplotlib==3.10.8 python3 -I -B -O"
$U code/check_passive_fermi.py > results/passive_fermi.csv
$U code/check_asymmetric_touch.py > results/asymmetric_touch.csv
$U code/certify_five_sites.py > results/five_site_certificate.txt
$U code/check_boundary_characters.py > results/ising_boundary.csv
$U code/reproduce_ising.py
$U code/ising_mesoscopic.py
$U code/crossover.py
$U code/dispersion.py
$U code/general_dispersion_check.py > results/general_dispersion.csv   # prints "all checks passed" on stderr
$U code/finite_size.py > results/finite_size.csv
python3 tools/compare_results.py <shipped results/> results/
```

Every script stops with "CHECK FAILED" and exit status 1 if a check fails; none uses `assert`, so they behave the
same under `python -O`. **Comparison policy:** `five_site_certificate.txt` (exact interval certificate) and
`ising_boundary.csv` (integer counts) must be byte-identical; in the other CSV files every non-numeric field must
match and every number must agree within relative 1e-6 (absolute 1e-7 for values below 1e-6, which are error
diagnostics or small differences), because floating-point and BLAS ordering change trailing digits across platforms
(on Linux CI the largest relative change in a reported value was about 2e-8, at L = 4096). Every number quoted in
the paper has at most seven significant digits. For finite_size.csv the tolerance is relative 1e-3: its L = 4096 entries for the cubic chain's one-body
control sum level energies of order 1e-10 from a dense eigensolver, and the column Q^G/Q changed from 1.09585 to
1.09568 between two Linux CI runners. The last digit of Table 4's entry 1.0959 is therefore platform-sensitive; the
value 1.096 quoted in the text, and all other Table 4 entries at their printed precision, are not. The figure PDFs embed
creation dates and are compared visually. On the release host every output except `ising_checks.csv` (relative
difference about 1e-15) was byte-identical to the shipped file. GitHub CI applies this policy on Linux and runs a
negative control.

## Licence

Text, paper, figures, results and records: CC0-1.0. Code: MIT. See [LICENSE.md](LICENSE.md).
