# Mathematics certificate: reproduction and release

This guide covers the computer-assisted claims in *The Separation of Peeling and Rank
Thresholds in Anisotropic Random Complexes*. It requires no external dataset.

## Clean environment and one-command driver

From the repository root:

```bash
python3 -m venv .venv-math
.venv-math/bin/python -m pip install -r requirements-math-lock.txt
.venv-math/bin/python scripts/reproduce_math_certificate.py
```

`requirements-math-lock.txt` pins `mpmath` and `python-flint`. The Bernstein coefficient and
rational Sturm checks use the Python standard library alone.

The driver runs the exact rational small-`r` lemma for `3 <= c <= 24`; the outward-rounded
interval certificate, which isolates both threshold roots for `3 <= c <= 9` and certifies
negative peeling surplus through `c = 24`; an independently implemented exact-polynomial
check using rational Sturm sequences for `3 <= c <= 9`; the exact Bernstein coefficient
certificate for `3 <= c <= 24`; and an independently implemented Arb spot check for
`c = 3, 4, 7, 24`. It writes:

- `reports/AUDIT_SMALL_R_LEMMA.md`;
- `reports/AUDIT_INTERVAL_CERTIFICATE.md`;
- `reports/AUDIT_POLYNOMIAL_CERTIFICATE.md`;
- `reports/AUDIT_ARB_SPOT_CHECK.md`;
- `reports/AUDIT_BERNSTEIN_CERTIFICATE.md`;
- `reports/MATH_CERTIFICATE_MANIFEST.json`.

The JSON manifest records Python and dependency versions, portable commands, scope verdicts,
SHA-256 hashes of every source file and certificate output, and one aggregate source-snapshot
digest. GitHub and the archival release record the root commit hash externally; a commit cannot
reliably contain its own hash.

The certificate scope is the declared code-only mathematics scope. Reproducing the certificate
depends on the six certificate scripts, the Markdown report helper, the dependency lock and
this guide. The manifest lists and hashes this complete scope and reports two verdicts:

- `scope_clean` — no uncommitted change in any path the certificate depends on.
- `worktree_clean` — no uncommitted change anywhere in the code-only repository.
- `releasable` — both verdicts are true.

A submission archive is valid when `releasable` is `true`. The scope is in the manifest so a
reader can verify exactly which sources determine the certificate. The generated manifest is
excluded from its own worktree verdict to avoid a circular self-check; the submission archive
and permanent release record its SHA-256 hash externally.

## Review and public release

This public repository contains the same code-only snapshot supplied as Supporting Information:
the scripts, dependency lock, generated reports, manifest and reproduction guide. Reviewers can
run either copy without access to the private development repository.

Before creating the submission archive:

1. commit every file in the declared certificate scope;
2. rerun the driver on that commit and verify `releasable: true`;
3. confirm that the manifest records the aggregate source digest and hashes every file in scope;
4. create the archive from that clean code snapshot.

The public repository contains certificate code, reports and reproducibility metadata. The
manuscript and the materials project remain outside it. Each published snapshot should carry a
version tag, the generated manifest and the full Git commit hash; a DOI archive may be added to
the release record when available.
