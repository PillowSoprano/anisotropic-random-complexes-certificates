"""Run every computer-assisted check used by the mathematics manuscript.

The release gate is the complete code-only certificate scope: six certificate scripts, one
report helper, the dependency lock and the reproduction guide. Files outside that declared
scope do not affect the executable certificate.

The manifest therefore records two verdicts. `scope_clean` is the release gate: no uncommitted
change in any path the certificate depends on. `worktree_clean` is reported alongside it for
information, and is deliberately not the gate. Every path in scope is listed in the manifest
and hashed, so a reader can check the scope rather than take it on trust.

Usage:  python3 scripts/reproduce_math_certificate.py
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "reports" / "MATH_CERTIFICATE_MANIFEST.json"
JOBS = [
    ("exact small-r lemma", ["scripts/run_small_r_lemma.py", "--c-max", "24"],
     "reports/AUDIT_SMALL_R_LEMMA.md", False),
    ("interval certificate", ["scripts/run_interval_certificate.py", "--c-max", "24",
                              "--threshold-c-max", "9"],
     "reports/AUDIT_INTERVAL_CERTIFICATE.md", False),
    ("independent polynomial check", ["scripts/run_polynomial_certificate.py"],
     "reports/AUDIT_POLYNOMIAL_CERTIFICATE.md", False),
    ("independent Arb spot checks", ["scripts/run_arb_spot_checks.py"],
     "reports/AUDIT_ARB_SPOT_CHECK.md", False),
    ("exact Bernstein coefficient check",
     ["scripts/verify_b7_sequences.py", "--c-max", "24"],
     "reports/AUDIT_BERNSTEIN_CERTIFICATE.md", True),
]

# Everything the certificate depends on. Directories are expanded to their tracked files.
# `src/prfield/report.py` is the Markdown report helper imported by the certificate scripts.
SCOPE = [
    "scripts/run_small_r_lemma.py",
    "scripts/run_interval_certificate.py",
    "scripts/run_polynomial_certificate.py",
    "scripts/run_arb_spot_checks.py",
    "scripts/verify_b7_sequences.py",
    "scripts/reproduce_math_certificate.py",
    "src/prfield/report.py",
    "requirements-math-lock.txt",
    "docs/MATH_CERTIFICATE_REPRODUCIBILITY.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git(*args: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True
    )
    return result.stdout.strip()


def main() -> None:
    try:
        dependency_versions = {
            "mpmath": importlib.metadata.version("mpmath"),
            "python-flint": importlib.metadata.version("python-flint"),
        }
    except importlib.metadata.PackageNotFoundError as error:
        raise SystemExit(
            "missing a locked certificate dependency; create the environment and install:\n"
            "  python3 -m venv .venv-math\n"
            "  .venv-math/bin/python -m pip install -r requirements-math-lock.txt"
        ) from error

    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    commands = []
    for label, arguments, report, capture in JOBS:
        command = [sys.executable, *arguments]
        portable_command = ["python", *arguments]
        print(f"\n[{label}] {' '.join(portable_command)}", flush=True)
        if capture:
            result = subprocess.run(
                command, cwd=ROOT, env=env, check=True, capture_output=True, text=True
            )
            (ROOT / report).write_text(
                "# Exact Bernstein coefficient certificate\n\n"
                "The table below is produced with exact rational arithmetic. A zero in the "
                "`eta recross` column and one in the `dq changes` column certify the "
                "single-crossing hypothesis used in Appendix B.\n\n```text\n"
                + result.stdout
                + "```\n"
            )
        else:
            subprocess.run(command, cwd=ROOT, env=env, check=True)
        report_path = ROOT / report
        if not report_path.is_file() or report_path.stat().st_size == 0:
            raise SystemExit(f"missing certificate output: {report}")
        commands.append({"label": label, "command": portable_command, "report": report})

    scope_files = sorted(git("ls-files", "--", *SCOPE).splitlines())
    if not scope_files:
        raise SystemExit("the declared scope matches no tracked file; refusing to certify")
    scope_dirty = sorted(git("status", "--porcelain", "--", *SCOPE).splitlines())
    files = sorted(set(scope_files) | {report for _, _, report, _ in JOBS})
    source_hashes = {name: sha256(ROOT / name) for name in scope_files}
    source_snapshot_sha256 = hashlib.sha256(
        json.dumps(source_hashes, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    record = {
        "scope": SCOPE,
        "scope_identity": "aggregate SHA-256 of the sorted source-file SHA-256 map",
        "source_snapshot_sha256": source_snapshot_sha256,
        "scope_clean": not scope_dirty,
        "scope_dirty_entries": scope_dirty,
        "worktree_clean": None,
        "dirty_entry_count": None,
        "releasable": None,
        "manifest_self_excluded_from_worktree_verdict": True,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": dependency_versions,
        "commands": commands,
        "sha256": {name: sha256(ROOT / name) for name in files},
    }
    # Write once, then audit every non-self-referential path. The second write records that
    # verdict; the archive supplies the manifest's own hash externally.
    MANIFEST.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    # The generated manifest cannot certify its own tracked bytes without a circular hash.
    # Audit every other path here; the archive/release supplies the manifest hash externally.
    dirty_paths = git(
        "status", "--porcelain", "--", ".", ":(exclude)reports/MATH_CERTIFICATE_MANIFEST.json"
    ).splitlines()
    record["worktree_clean"] = not dirty_paths
    record["dirty_entry_count"] = len(dirty_paths)
    record["releasable"] = not scope_dirty and not dirty_paths
    MANIFEST.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {MANIFEST.relative_to(ROOT)}")
    if scope_dirty:
        print("NOT RELEASABLE: uncommitted changes inside the certificate scope:")
        for line in scope_dirty:
            print(f"  {line}")
        print("Commit them and rerun before archiving.")
    else:
        print(f"scope_clean: true over {len(scope_files)} tracked files")
    if dirty_paths:
        print(f"NOT RELEASABLE: {len(dirty_paths)} uncommitted entries remain in the worktree")
    else:
        print("worktree_clean: true; releasable: true")


if __name__ == "__main__":
    main()
