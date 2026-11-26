# Replay receipt (v0.1.0-candidate, 26 November 2026)

Release replay by the publisher from a fresh copy of the approved candidate's code/, macOS arm64, working directory =
package root, one CPU, about one minute in all, OPENBLAS_NUM_THREADS=1, MPLBACKEND=Agg. Interpreter as printed inside
the run (see [ENVIRONMENT.txt](ENVIRONMENT.txt)): CPython 3.12.11 with numpy 2.3.5, scipy 1.17.0, mpmath 1.3.0,
matplotlib 3.10.8.

| Program | Flags | Exit | Output |
|---|---|---|---|
| check_passive_fermi.py | -I -B -O | 0 | byte-identical to results/passive_fermi.csv |
| check_asymmetric_touch.py | -I -B -O | 0 | byte-identical to results/asymmetric_touch.csv |
| certify_five_sites.py | -I -B -O | 0 | byte-identical to results/five_site_certificate.txt |
| check_boundary_characters.py | -I -B -O | 0 | byte-identical to results/ising_boundary.csv |
| reproduce_ising.py | -I -B -O | 0 | results/ising_checks.csv within the comparison policy (final digits only) |
| ising_mesoscopic.py | -I -B -O | 0 | byte-identical to results/ising_mesoscopic.csv |
| crossover.py | -I -B -O | 0 | byte-identical to results/crossover.csv; figures/crossover.pdf regenerated |
| dispersion.py | -I -B -O | 0 | figures/dispersion.pdf regenerated |
| general_dispersion_check.py | -I -B -O | 0 | byte-identical to results/general_dispersion.csv; stderr "all checks passed" |
| finite_size.py | -I -B -O | 0 | byte-identical to results/finite_size.csv |
| Negative control | general_dispersion_check.py with expected D_3 changed to 11.3023136908, -I -B -O | 1 | "CHECK FAILED: D_3" |

`tools/compare_results.py` passed on all nine result files. The outputs match the lead's replay of 25 November 2026.
After the replay only the title-page DOI line of manuscript.tex, the rebuilt PDF and the release records changed;
code/, results/ and figures/ are byte-identical to the approved candidate. GitHub CI on Linux reruns all scripts under
the comparison policy in ENVIRONMENT.txt, with the negative control.
