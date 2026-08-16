# Certificates for anisotropic random complexes

This code-only repository contains the exact and interval-arithmetic certificates for
*The Separation of Peeling and Rank Thresholds in Anisotropic Random Complexes*. It contains
no manuscript, proof appendix, PDF, materials-project file or submission source.

## Certified statements

- Exact Bernstein coefficients certify the shape condition for `3 <= c <= 24`.
- Outward-rounded interval arithmetic certifies strict peeling/Maxwell separation for
  `3 <= c <= 24`.
- For `3 <= c <= 9`, the interval calculation isolates both threshold roots.
- An independent rational Sturm calculation verifies the shape condition for `3 <= c <= 9`.
- A separately implemented `python-flint`/Arb calculation checks the complete peeling boxes
  for `c = 3, 4, 7, 24` and the tabulated Maxwell boxes for `c = 3, 4, 7`.

The uniform all-`c` shape condition remains open.

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
