"""The last gap: certify D_c < 0 on (0, r0], where no covering can work.

`D_c(0) = 0` to second order, so interval covering cannot certify a sign in a right
neighbourhood of the origin. That stretch needs a series argument, and this is it.

`H_c` is a polynomial with exact rational coefficients (from the connectivity recursion), and
`s(r)(1-r) = r - sum_{j>=2} r^j / (j(j-1))` has coefficients bounded by `1/2` in modulus. So
`D_c = H_c - s(1-r)H_c'` has computable rational Taylor coefficients and an explicit tail
bound: truncating at order `N`,

    |tail| <= (||H_c'||_1 / 2) * r^{N+1} / (1 - r),

whence on `(0, r0]`

    D_c(r) / r^2  <=  d_2 + sum_{k=3..N} |d_k| r0^{k-2} + (||H_c'||_1 / 2) r0^{N-1} / (1 - r0),

and the right-hand side is a single exact rational. Negative means `D_c < 0` on the whole of
`(0, r0]`.

The leading coefficient comes out as `d_2 = -(2c - 5)/2`, exactly, so it is at most `-1/2` for
every `c >= 3` -- the sign of the whole argument is structural, not numerical.

Usage:  PYTHONPATH=src python3 scripts/run_small_r_lemma.py --c-max 24
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from run_interval_certificate import H_poly, p_deriv, p_mul  # noqa: E402

from prfield.report import write_report  # noqa: E402

REPORT = Path("reports/AUDIT_SMALL_R_LEMMA.md")
N = 48
R0 = F(1, 200)


def s_times_one_minus_r(n: int) -> list[F]:
    """Taylor coefficients of s(r)(1-r) with s = -log(1-r), to order n."""
    c = [F(0)] * (n + 1)
    if n >= 1:
        c[1] = F(1)
    for j in range(2, n + 1):
        c[j] = -F(1, j * (j - 1))
    return c


def analyse(c: int, n: int = N, r0: F = R0) -> dict:
    Hp = H_poly(c)
    Hd = p_deriv(Hp)
    prod = p_mul(s_times_one_minus_r(n), Hd)
    Hx = (Hp + [F(0)] * (n + 1))[: n + 1]
    px = (prod + [F(0)] * (n + 1))[: n + 1]
    D = [Hx[k] - px[k] for k in range(n + 1)]
    norm = sum(abs(x) for x in Hd)
    bound = (D[2]
             + sum(abs(D[k]) * r0 ** (k - 2) for k in range(3, n + 1))
             + F(norm, 2) * r0 ** (n - 1) / (1 - r0))
    return {"c": c, "d0": D[0], "d1": D[1], "d2": D[2],
            "predicted_d2": -F(2 * c - 5, 2), "norm": norm,
            "bound": bound, "certified": bound < 0}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-min", type=int, default=3)
    parser.add_argument("--c-max", type=int, default=9)
    parser.add_argument("--report", type=Path, default=REPORT)
    args = parser.parse_args()
    if args.c_min < 3 or args.c_max < args.c_min:
        raise SystemExit("require 3 <= --c-min <= --c-max")

    rows = [analyse(c) for c in range(args.c_min, args.c_max + 1)]
    lines = ["# The small-`r` lemma\n"]
    lines.append(
        "`D_c(0) = 0` to second order, so no interval covering can certify a sign near the "
        "origin. This closes `(0, r0]` by an exact series argument instead, with "
        f"`r0 = {float(R0)}` and truncation order `N = {N}`.\n"
    )
    lines.append("| c | d₀ | d₁ | d₂ (exact) | −(2c−5)/2 | ‖H_c'‖₁ | bound on D_c/r² | certified |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for r in rows:
        lines.append(
            f"| {r['c']} | {r['d0']} | {r['d1']} | {r['d2']} | {r['predicted_d2']} | "
            f"{float(r['norm']):.3g} | {float(r['bound']):.6f} | "
            f"{'yes' if r['certified'] else 'NO'} |")
    lines.append(
        "\n## What this establishes\n"
        "`d₀ = d₁ = 0` and `d₂ = −(2c−5)/2` exactly, so the vanishing at the origin is exactly "
        "second order and the leading coefficient is at most `−1/2` for every `c ≥ 3`. The "
        "bound column is a single exact rational upper bound for `D_c(r)/r²` on the whole of "
        "`(0, r0]`, computed from the exact Taylor coefficients plus the explicit tail bound; "
        "being negative, it certifies `D_c < 0` there.\n"
        "\nFor the globally covered cases this closes the interval certificate at the origin. "
        "Across the full tested range it independently confirms the left-endpoint sign used by "
        "the exact Bernstein single-crossing certificate.\n"
        "\n## What it does not\n"
        f"This is a per-instance certificate for `{args.c_min} ≤ c ≤ {args.c_max}`. A uniform "
        "proof for all `c ≥ 3` requires a bound on the Taylor coefficients that is independent "
        "of `c`.\n"
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_report(args.report, "\n".join(lines))
    for r in rows:
        print(f"c={r['c']:>2}  d2={str(r['d2']):>6} (predicted {r['predicted_d2']})  "
              f"bound={float(r['bound']):+.6f}  {'certified' if r['certified'] else 'FAILED'}")
    print(f"\nwrote {args.report}")
    if not all(r["certified"] for r in rows):
        raise SystemExit("small-r certificate failed for one or more requested c values")


if __name__ == "__main__":
    main()
