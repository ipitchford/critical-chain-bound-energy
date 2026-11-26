# Review summary (internal, AI-assisted)

This is a public summary of the internal checking of version 0.1.0-candidate, without scores. It is not peer review,
independent reproduction or formal verification. No human referee has assessed the paper.

## What was done (24-25 November 2026)

1. **A Claude review panel** in four roles (editor, domain specialist, perspective, devil's advocate). No reviewer
   found a mathematical error. The panel asked for clearer framing of what the conjecture says, accurate credit to
   earlier work, finite-size evidence for the limits, and a more readable structure.
2. **A methodology review with a line-by-line proof audit** and independent code. No error in any proof; minor gaps
   in written steps and one framing point (limits must not be read as finite-size ratios).
3. **A clean replay** of all scripts from a fresh copy; reported values reproduced.
4. **One revision round with extensions**: the three forms of the conjecture separated, with a status table; the
   rank-matching lemma stated in general and used by every proof; a general theorem for every sign-preserving
   polynomial dispersion, with the "if and only if" condition for the conformal law; a finite-size table up to
   L = 4096 with its script; new references (Allahverdyan-Balian-Nieuwenhuizen, Calabrese-Lefevre,
   Deift-Its-Krasovsky, Goldstein-Sela), all obtained and read; the one external input named as such; all scripts
   converted to explicit checks that also fail under `python -O`.
5. **Re-review and proof audit of the new and rewritten proofs.** No error. Remaining items (velocity normalisation in
   the comparison with Mitra and Srivastava, wording of the charge-sector floor, two numbers corrected against fresh
   outputs) were applied before release.

No Astra/Codex second-model review took place in this batch. Attempted but not landed, and stated as open problems:
a test on the interacting Heisenberg/XXZ chain, a lattice termination with a fixed boundary condition, next-order
finite-size corrections.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| The conjecture was treated as one statement | Three forms separated: (A) all critical states, false; (B) CFT-described states, proved for regulated CFT spectra and the hopping and Ising chains, open for interacting chains; (C) slope 6/π, true for unit hopping only (Table 5) |
| Mula et al.'s formula is a one-body optimum, not the unrestricted one | Stated in Sections 1-2; Theorem 11 at z = 1 shows that both readings agree for hopping |
| Credit to Mitra and Srivastava incomplete | Their analytic and numerical findings described as published, including the deviation at all-to-all pairing |
| The lattice proofs rest on an outside result | The boundary Rényi asymptotics of the open hopping chain named as the one external input (abstract, Appendix A, Discussion) |
| Limits presented without finite-size evidence | Table 4 and `code/finite_size.py`; every coefficient labelled as an L → ∞ limit; slow convergence stated |
| Ising free/fixed example overstated | Labelled as a continuum illustration with its assumption stated |
| Scripts used `assert` | Explicit checks with exit status 1, active under `-O` |

## Limits that remain

Free-fermion chains and an idealised conformal model only; interacting lattice chains remain open. One external
input (boundary Rényi asymptotics of the open hopping chain). Idealised operations without limits on depth or
locality. Coefficients are limits as L grows; some converge only logarithmically slowly. The computations check the
stated finite-size values, constants and identities; the theorems rest on the written proofs.
