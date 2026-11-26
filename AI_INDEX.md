# AI index: The energy an isolated half of a critical chain cannot give up (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1
states the problem, why it matters, what was known, what is new, what is not claimed and a roadmap; Section 9 gives
the status of the conjecture and open problems; the "Data and code availability" section describes the scripts).
DOI 10.5281/zenodo.24553557 (concept 10.5281/zenodo.24553556).
Status: unrefereed candidate; internal AI-assisted review, proof audits and computation only.

## Notation

Ground state of a chain of N = 2L sites; ρ_L reduced state of the left half; H_L block Hamiltonian (terms crossing the
cut removed). Q_L = bound (residual) energy: min over unitaries U of Tr(U ρ_L U* H_L) − E_0(H_L), i.e. the passive
energy above the block ground energy. Q_L^ch: minimum over particle-number-conserving unitaries; Q_L^G: minimum over
one-body (Gaussian) number-conserving rotations. S_L = entanglement entropy of ρ_L. v velocity, c central charge.
T_{2L} open hopping matrix; dΓ(f(T)) free-fermion chain with one-particle Hamiltonian f(T); m = largest order of a
zero of f on [−2, 2]; d = log L.

## Claims and where they are proved

| ID | Claim | Location | Status / dependence |
|---|---|---|---|
| RM | Rank matching: under hypotheses (a)-(c), Q_L/g_L → φ(1) (isentropic reading: energy of the block thermal state with entropy S_L) | Lemma 1, Remark 2 | written proof |
| CFT | Regulated boundary CFT: Q_L ~ (π v/(4 K_p L)) S_L²; vacuum-dominated unitary CFT: Q_L ~ 3v S_L²/(2π c L) | Theorem 3 | written proof (model of the reduced state) |
| BC | Next term fixed by the boundary entropy of the physical termination (rational CFT); passive state → strip thermal state in relative entropy | Theorem 4, Corollary 5 | written proof |
| GEN | Every sign-preserving polynomial f: dΓ(f(T_{2L})) has the hopping ground state, S_L = (1/6) log L + O(1), Q_L ~ D(m,A) S_L^{m+1}/L^m | Theorem 6 | written proof; uses the boundary Rényi asymptotics of the open hopping chain (Appendix A; Fagotti-Calabrese Eq. (38), from Deift-Its-Krasovsky) |
| ODD | f(x) = x^z, odd z: D_1 = 3/π, D_3 = 11.2023137...; LQ_L/S_L² → 0 for z > 1 | Corollary 7 | written proof; constants checked by general_dispersion_check.py |
| IFF | L Q_L/S_L² has a positive finite limit iff f'(0) ≠ 0 and f has no other zero in [−2, 2]; then Q_L ~ 3v S_L²/(2π L), v = 2|f'(0)| | Corollary 8 | written proof |
| ISING | Critical Ising chain: S_I = (1/12) log L + O(1), Q_I ~ (6/π) S_I²/L, N Q_I/S_I² → 12/π | Corollary 9 | written proof; spin vs Majorana checks |
| CROSS | Linear perturbation of a flat dispersion: explicit crossover function e_z(λ), asymptote restoring S²/L | Theorem 10 | written proof; quadrature in crossover.py |
| GAUSS | One-body control for f = x^z: Q_L^G limit ratio R_z to Q_L (R_1 = 1, R_3 ≈ 2.14815, R_5 ≈ 11.55617) | Theorem 11 | written proof; constants checked |
| CH | Charge conservation costs nothing at leading order for f = x^z | Corollary 12, Appendix C | written proof |
| THREE | Range-three chain f(x) = x(x − 3/2)²: Q_L, Q_L^ch, Q_L^G of orders d³/L², d/L, d²/L | Theorem 13, Table 2 | written proof |
| PQ | Unequal particle and hole orders f = x^p (x − a)^q | Corollary 14 | written proof |

Conjecture status (Table 5): form (A) "all critical states" false; form (B) "CFT-described states" proved for
regulated CFT spectra and the open hopping and Ising chains, open for interacting lattice chains; form (C) slope 6/π
true for unit hopping only; general slope 3v/(π c) in N units.

## Evidence

| Check | Program | Output | Arithmetic |
|---|---|---|---|
| Complete Fock-space ordering and direct many-body checks, odd powers (L = 6, 16) | [code/check_passive_fermi.py](code/check_passive_fermi.py) | [results/passive_fermi.csv](results/passive_fermi.csv) | floats |
| Range-three chain: Q, Q^ch, Q^G at small L | [code/check_asymmetric_touch.py](code/check_asymmetric_touch.py) | [results/asymmetric_touch.csv](results/asymmetric_touch.csv) | floats |
| Ordering inversion at N = 10, L = 5 (ratio about 34.84, gap 2 − √3) | [code/certify_five_sites.py](code/certify_five_sites.py) | [results/five_site_certificate.txt](results/five_site_certificate.txt) | mpmath interval arithmetic |
| Ising annulus spectra, exact multiplicities (Table 3 right) | [code/check_boundary_characters.py](code/check_boundary_characters.py) | [results/ising_boundary.csv](results/ising_boundary.csv) | integers and floats |
| Ising spin vs Majorana agreement (within 2e-14) | [code/reproduce_ising.py](code/reproduce_ising.py) | [results/ising_checks.csv](results/ising_checks.csv) | floats |
| Ising mesoscopic values (Table 3 left) | [code/ising_mesoscopic.py](code/ising_mesoscopic.py) | [results/ising_mesoscopic.csv](results/ising_mesoscopic.csv) | floats |
| Crossover e_3(λ), Figure 1 | [code/crossover.py](code/crossover.py) | [results/crossover.csv](results/crossover.csv), [figures/crossover.pdf](figures/crossover.pdf) | quadrature |
| Dispersion figure (Figure 2) | [code/dispersion.py](code/dispersion.py) | [figures/dispersion.pdf](figures/dispersion.pdf) | plot |
| Pressure constant A (four polynomials), D_1, D_3, range-three coefficient 2.0868, R_3, R_5, large-z ratio, pair ratio 1.25734 | [code/general_dispersion_check.py](code/general_dispersion_check.py) | [results/general_dispersion.csv](results/general_dispersion.csv) | mpmath, floats |
| Finite-size table L = 16 ... 4096 (Table 4) | [code/finite_size.py](code/finite_size.py) | [results/finite_size.csv](results/finite_size.csv) | floats |

Negative control: changing the expected D_3 in general_dispersion_check.py makes it exit 1 under python -O with
"CHECK FAILED: D_3" ([REPLAY_RECEIPT.md](REPLAY_RECEIPT.md); also run in CI). The computations check finite-size
values, constants and identities; they do not certify the asymptotic theorems, which rest on the written proofs.

Numbers quoted on the release page: 3v/(2πc); slopes 3/π (hopping) and 6/π (Ising) in L units, 6/π and 12/π in N
units; charge-conserving floor over Q_L about 20 times at L = 16 and about 1,360 times at L = 4096; cubic chain
L³Q_L/S_L⁴ = 36.3 at L = 4096 against the limit 11.2, and Q_L^G/Q_L = 1.096 at L = 4096 against the limit 2.148;
orders (log L)³/L², (log L)/L, (log L)²/L.

## Limits

Free-fermion lattice chains and a regulated conformal model only; interacting chains (e.g. Heisenberg) open. One
external input (boundary Rényi asymptotics of the open hopping chain). Idealised unitaries (no depth or locality
limit). All coefficients are L → ∞ limits; some converge logarithmically slowly (Table 4). Literature search of
limited scope.

## Replay and records

Replay commands: [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml), comparator [tools/compare_results.py](tools/compare_results.py).
Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json).
Bibliography: [refs.bib](refs.bib). Licences: [LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes:
[MANIFEST.sha256](MANIFEST.sha256).
