# Certificates for anisotropic random complexes

**Current manuscript (30 September 2026): _Exact Rank and Maxwell Thresholds in Anisotropic Random Complexes_.**
This repository is the public base snapshot. The manuscript's Supporting Information
(`AAP_Supporting_Information.zip`, indexed by `CERTIFICATE_INDEX.json`) includes the
extended shape certificates for `3 <= c <= 200` and the global-maximiser certificates
that complete the rank result for every fixed `c >= 3`. The uniform all-`c` shape
condition remains open (see Conjecture 80); the all-`c` rank formula (Theorem 1) is
unconditional, using the global-maximiser selection of Theorem 51.

This code-only base snapshot was originally released for
*The Separation of Peeling and Rank Thresholds in Anisotropic Random Complexes*. It contains
no manuscript, proof appendix, PDF, materials-project file or submission source.

## Certified statements in this base snapshot

- Exact Bernstein coefficients certify the shape condition for `3 <= c <= 24`.
- Outward-rounded interval arithmetic certifies strict peeling/Maxwell separation for
  `3 <= c <= 24`.
- For `3 <= c <= 9`, the interval calculation isolates both threshold roots.
- An independent rational Sturm calculation verifies the shape condition for `3 <= c <= 9`.
- A separately implemented `python-flint`/Arb calculation checks the complete peeling boxes
  for `c = 3, 4, 7, 24` and the tabulated Maxwell boxes for `c = 3, 4, 7`.

## Reproduce

```bash
python3 -m venv .venv-math
.venv-math/bin/python -m pip install -r requirements-math-lock.txt
.venv-math/bin/python scripts/reproduce_math_certificate.py
```

The driver writes five reports and `reports/MATH_CERTIFICATE_MANIFEST.json`. The manifest
records the commands, versions, clean-scope verdicts, file hashes and aggregate source-snapshot
digest.

See `docs/MATH_CERTIFICATE_REPRODUCIBILITY.md` for the complete protocol. This public
repository was created from the audited code-only snapshot supplied with the manuscript. The
private development repository, manuscript source and materials project are not part of this
release.

## Citation and licence

The repository is released under the MIT License. Cite the archived release metadata in
`CITATION.cff` and record the full Git commit when reporting a reproduction run.
