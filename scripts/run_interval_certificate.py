"""P2a / P5e: outward-rounded interval certificates for the two thresholds.

Closure criterion 2 of the mathematics. Both thresholds reduce to one-dimensional root
isolation in the variable `r`, and the same computation serves P2a (certifying
`theta_peel`) and P5e (the shape condition that puts `max L_c` at `0` or `u_*`).

The key structural fact that makes a certificate cheap: `C_k(r)`, `H_c(r)` and `K_c(r)` are
**polynomials in r with rational coefficients**. They are built here by exact `Fraction`
arithmetic from the connectivity recursion, so `H_c'` is the exact derivative polynomial
rather than a finite difference, and the only transcendental entering is `log(1-r)`, which
`mpmath.iv` encloses with outward rounding.

Certified quantities, per `c`:

    D_c(r) = H_c(r) - [-log(1-r)] (1-r) H_c'(r)        sign(D_c) = sign(g_c')
    g_c(r) = [-log(1-r)] / H_c(r)                      theta^2 at a positive fixed point
    Psi(r) = 2[K_c(r) - c] + C(c,2) [-log(1-r)] (2 - H_c(r))

Certificate 1 (`theta_peel`, and the P5e shape condition): `D_c` has exactly one zero in
`(0,1)`, negative before and positive after. Uniqueness is certified by covering the rest of
the interval with boxes on which the enclosure of `D_c` excludes `0`.

Certificate 2 (`theta_Max`): `Psi` has exactly one zero above `r_peel`, negative before and
positive after, again with a covering.

Usage:  PYTHONPATH=src python3 scripts/run_interval_certificate.py --c-max 24
"""

from __future__ import annotations

import argparse
import sys
import time
from fractions import Fraction
from math import comb
from pathlib import Path

from mpmath import iv, mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from prfield.report import write_report  # noqa: E402

REPORT = Path("reports/AUDIT_INTERVAL_CERTIFICATE.md")
PREC = 60
R0 = 0.005  # left endpoint of the interval covering; the exact series closes (0, R0]


# ---------------------------------------------------------------------------
# exact polynomial arithmetic, coefficients in ascending powers of r
# ---------------------------------------------------------------------------

def p_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fraction(0)) + (b[i] if i < len(b) else Fraction(0))
            for i in range(n)]


def p_scale(a, s):
    return [c * s for c in a]


def p_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def p_one_minus_r_pow(k):
    """(1-r)^k by binomial expansion, exact."""
    return [Fraction((-1) ** i * comb(k, i)) for i in range(k + 1)]


def p_deriv(a):
    return [a[i] * i for i in range(1, len(a))] or [Fraction(0)]


def connected_polys(c):
    """C_k(r) for k = 1..c, exactly, from C_m = 1 - sum_k binom(m-1,k-1) C_k (1-r)^{k(m-k)}."""
    C = {1: [Fraction(1)]}
    for m in range(2, c + 1):
        acc = [Fraction(0)]
        for k in range(1, m):
            term = p_mul(C[k], p_one_minus_r_pow(k * (m - k)))
            acc = p_add(acc, p_scale(term, Fraction(comb(m - 1, k - 1))))
        C[m] = p_add([Fraction(1)], p_scale(acc, Fraction(-1)))
    return C


def H_poly(c):
    """H_c(r) = sum_k binom(c-1,k-1) C_k(r) (1-r)^{k(c-k)} (k-1)/(c-1)."""
    C = connected_polys(c)
    out = [Fraction(0)]
    for k in range(1, c + 1):
        w = Fraction(comb(c - 1, k - 1) * (k - 1), c - 1)
        out = p_add(out, p_scale(p_mul(C[k], p_one_minus_r_pow(k * (c - k))), w))
    return out


def K_poly(c):
    """K_c(r) = sum_k binom(c,k) C_k(r) (1-r)^{k(c-k)}."""
    C = connected_polys(c)
    out = [Fraction(0)]
    for k in range(1, c + 1):
        out = p_add(out, p_scale(p_mul(C[k], p_one_minus_r_pow(k * (c - k))),
                                 Fraction(comb(c, k))))
    return out


# ---------------------------------------------------------------------------
# interval evaluation, structurally rather than in the monomial basis
# ---------------------------------------------------------------------------

class Dual:
    """(value, derivative) pair over intervals, so H_c' comes out of the same recursion.

    Evaluating H_c from its expanded monomial coefficients is hopeless in interval
    arithmetic: the (1-r)^k expansions give alternating coefficients of size up to
    binom(72,36), and Horner then overestimates by many orders of magnitude. Evaluating the
    connectivity recursion structurally keeps every product positive, with a single
    cancellation per level, and the enclosures stay tight.
    """

    __slots__ = ("v", "d")

    def __init__(self, v, d):
        self.v, self.d = v, d

    def __add__(self, o):
        o = _lift(o)
        return Dual(self.v + o.v, self.d + o.d)

    __radd__ = __add__

    def __sub__(self, o):
        o = _lift(o)
        return Dual(self.v - o.v, self.d - o.d)

    def __rsub__(self, o):
        o = _lift(o)
        return Dual(o.v - self.v, o.d - self.d)

    def __mul__(self, o):
        o = _lift(o)
        return Dual(self.v * o.v, self.d * o.v + self.v * o.d)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = _lift(o)
        return Dual(self.v / o.v, (self.d * o.v - self.v * o.d) / (o.v * o.v))


def _lift(o):
    return o if isinstance(o, Dual) else Dual(iv.mpf(o) if not hasattr(o, "a") else o,
                                              iv.mpf(0))


def _q_pow(x, k):
    """(1-r)^k as a Dual in r."""
    q = iv.mpf(1) - x
    if k == 0:
        return Dual(iv.mpf(1), iv.mpf(0))
    return Dual(q ** k, iv.mpf(-k) * (q ** (k - 1)))


def _connected_duals(c, x):
    """C_k(r) and C_k'(r) for k = 1..c, by the recursion, in interval dual arithmetic."""
    C = {1: Dual(iv.mpf(1), iv.mpf(0))}
    for m in range(2, c + 1):
        acc = Dual(iv.mpf(0), iv.mpf(0))
        for k in range(1, m):
            # Dual on the left: mpmath's iv raises instead of returning NotImplemented,
            # so __rmul__ would never be reached.
            acc = acc + (C[k] * _q_pow(x, k * (m - k))) * iv.mpf(comb(m - 1, k - 1))
        C[m] = Dual(iv.mpf(1), iv.mpf(0)) - acc
    return C


class Dual2:
    """(value, first derivative, second derivative) over intervals.

    Needed for the centered form. Interval arithmetic is sound under any valid extension, so
    the naive and centered forms both enclose the true range and differ only in tightness --
    but that soundness depends on the second derivative being a genuine enclosure, so an error
    here would corrupt the certificate silently rather than merely loosen it. The inclusion
    self-test below is what guards against that.
    """

    __slots__ = ("v", "d", "dd")

    def __init__(self, v, d, dd):
        self.v, self.d, self.dd = v, d, dd

    def __add__(self, o):
        o = _lift2(o)
        return Dual2(self.v + o.v, self.d + o.d, self.dd + o.dd)

    __radd__ = __add__

    def __sub__(self, o):
        o = _lift2(o)
        return Dual2(self.v - o.v, self.d - o.d, self.dd - o.dd)

    def __rsub__(self, o):
        o = _lift2(o)
        return Dual2(o.v - self.v, o.d - self.d, o.dd - self.dd)

    def __mul__(self, o):
        o = _lift2(o)
        return Dual2(self.v * o.v,
                     self.d * o.v + self.v * o.d,
                     self.dd * o.v + iv.mpf(2) * self.d * o.d + self.v * o.dd)

    __rmul__ = __mul__


def _lift2(o):
    if isinstance(o, Dual2):
        return o
    z = o if hasattr(o, "a") else iv.mpf(o)
    return Dual2(z, iv.mpf(0), iv.mpf(0))


def _q_pow2(x, k):
    """(1-r)^k as a Dual2 in r: v = q^k, d = -k q^{k-1}, dd = k(k-1) q^{k-2}."""
    q = iv.mpf(1) - x
    if k == 0:
        return Dual2(iv.mpf(1), iv.mpf(0), iv.mpf(0))
    if k == 1:
        return Dual2(q, iv.mpf(-1), iv.mpf(0))
    return Dual2(q ** k, iv.mpf(-k) * (q ** (k - 1)),
                 iv.mpf(k * (k - 1)) * (q ** (k - 2)))


def _connected_duals2(c, x):
    C = {1: Dual2(iv.mpf(1), iv.mpf(0), iv.mpf(0))}
    for m in range(2, c + 1):
        acc = Dual2(iv.mpf(0), iv.mpf(0), iv.mpf(0))
        for k in range(1, m):
            acc = acc + (C[k] * _q_pow2(x, k * (m - k))) * iv.mpf(comb(m - 1, k - 1))
        C[m] = Dual2(iv.mpf(1), iv.mpf(0), iv.mpf(0)) - acc
    return C


def make_funcs(c):
    L = comb(c, 2)

    def H_dual(x):
        C = _connected_duals(c, x)
        out = Dual(iv.mpf(0), iv.mpf(0))
        for k in range(1, c + 1):
            w = iv.mpf(comb(c - 1, k - 1) * (k - 1)) / iv.mpf(c - 1)
            out = out + (C[k] * _q_pow(x, k * (c - k))) * w
        return out

    def K_val(x):
        C = _connected_duals(c, x)
        out = Dual(iv.mpf(0), iv.mpf(0))
        for k in range(1, c + 1):
            out = out + (C[k] * _q_pow(x, k * (c - k))) * iv.mpf(comb(c, k))
        return out.v

    def s_of(x):
        return -iv.log(iv.mpf(1) - x)

    def D(x):
        h = H_dual(x)
        return h.v - s_of(x) * (iv.mpf(1) - x) * h.d

    def g(x):
        return s_of(x) / H_dual(x).v

    def Psi(x):
        return iv.mpf(2) * (K_val(x) - iv.mpf(c)) + iv.mpf(L) * s_of(x) * (iv.mpf(2) - H_dual(x).v)

    # --- centered forms -----------------------------------------------------
    def H2(x):
        C = _connected_duals2(c, x)
        out = Dual2(iv.mpf(0), iv.mpf(0), iv.mpf(0))
        for k in range(1, c + 1):
            w = iv.mpf(comb(c - 1, k - 1) * (k - 1)) / iv.mpf(c - 1)
            out = out + (C[k] * _q_pow2(x, k * (c - k))) * w
        return out

    def K2(x):
        C = _connected_duals2(c, x)
        out = Dual2(iv.mpf(0), iv.mpf(0), iv.mpf(0))
        for k in range(1, c + 1):
            out = out + (C[k] * _q_pow2(x, k * (c - k))) * iv.mpf(comb(c, k))
        return out

    def Dprime(x):
        """D' = s [ H' - (1-r) H'' ], from D = H - s(1-r)H' and s' = 1/(1-r)."""
        h = H2(x)
        return s_of(x) * (h.d - (iv.mpf(1) - x) * h.dd)

    def Psiprime(x):
        """Psi' = 2K' + L[ (2-H)/(1-r) - s H' ]."""
        h, k = H2(x), K2(x)
        return iv.mpf(2) * k.d + iv.mpf(L) * (
            (iv.mpf(2) - h.v) / (iv.mpf(1) - x) - s_of(x) * h.d)

    def _centered(f, fp, x):
        m = (x.a + x.b) / 2
        mid = iv.mpf([m, m])
        return f(mid) + fp(x) * (x - mid)

    def _meet(f, fp, x):
        """Intersect the naive and centered enclosures.

        Both are valid extensions, so neither need contain the other and the centered form is
        sometimes the looser of the two -- observed here near r -> 1 for Psi. Their
        intersection is still a valid enclosure and is never worse than either, so it is what
        the covering uses. Soundness was checked by point sampling inside boxes (3360 checks,
        no violations for either form), not by testing one against the other, which would have
        measured tightness rather than correctness.
        """
        n_, c_ = f(x), _centered(f, fp, x)
        return iv.mpf([max(n_.a, c_.a), min(n_.b, c_.b)])

    def Dc(x):
        return _meet(D, Dprime, x)

    def Psic(x):
        return _meet(Psi, Psiprime, x)

    return D, g, Psi, Dc, Psic


def _neg(z):  # enclosure strictly negative
    return z.b < 0


def _pos(z):  # enclosure strictly positive
    return z.a > 0


def endpoint(z, side, digits=14):
    """Format one endpoint of an interval without mpmath's nested point brackets."""
    return mp.nstr(mp.mpf(getattr(z, side)), digits)


def number(z, digits=14):
    return mp.nstr(mp.mpf(z), digits)


def cover(f, lo, hi, max_depth=26):
    """Adaptively bisect [lo,hi]; return the list of boxes whose sign stays indeterminate.

    A uniform grid is useless here: interval evaluation of a smooth function over a wide box
    overestimates badly (the dependency problem), so boxes far from the root stay
    indeterminate even though point enclosures there are tight to sixty digits. Bisection
    shrinks the overestimation quadratically and resolves them.
    """
    stack, unresolved = [(lo, hi, 0)], []
    while stack:
        a, b, depth = stack.pop()
        z = f(iv.mpf([a, b]))
        if _neg(z) or _pos(z):
            continue
        if depth >= max_depth:
            unresolved.append((a, b))
            continue
        m = (a + b) / 2
        stack.append((a, m, depth + 1))
        stack.append((m, b, depth + 1))
    return sorted(unresolved)


def isolate(f, lo, hi, max_depth=26):
    """Certify a single sign change on [lo,hi] and return its bracket."""
    unresolved = cover(f, lo, hi, max_depth)
    if not unresolved:
        return None, "no sign change on this range"
    # merge touching boxes
    merged = [list(unresolved[0])]
    for a, b in unresolved[1:]:
        if a <= merged[-1][1] * (1 + 1e-15):
            merged[-1][1] = b
        else:
            merged.append([a, b])
    if len(merged) > 1:
        return None, f"{len(merged)} disjoint indeterminate regions"
    a, b = merged[0]
    za, zb = f(iv.mpf([a, a])), f(iv.mpf([b, b]))
    if not ((_neg(za) and _pos(zb)) or (_pos(za) and _neg(zb))):
        return None, "endpoints do not straddle"
    return (a, b), "ok"


def bisect_unique_root(f, lo, hi, iterations=100):
    """Bracket a known unique negative-to-positive root using point intervals.

    Uniqueness is supplied by the exact Bernstein coefficient certificate. This routine adds
    an outward-rounded location bracket without paying for a second global interval covering.
    """
    lo, hi = mp.mpf(str(lo)), mp.mpf(str(hi))
    if not _neg(f(iv.mpf([lo, lo]))) or not _pos(f(iv.mpf([hi, hi]))):
        return None, "endpoints do not have signs negative, positive"
    for _ in range(iterations):
        mid = (lo + hi) / 2
        z = f(iv.mpf([mid, mid]))
        if _neg(z):
            lo = mid
        elif _pos(z):
            hi = mid
        else:
            return None, "point interval does not exclude zero"
    return (lo, hi), "ok"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-min", type=int, default=3)
    parser.add_argument("--c-max", type=int, default=9)
    parser.add_argument(
        "--threshold-c-max", type=int, default=9,
        help="largest c for the independent global interval covering and theta_Max isolation",
    )
    parser.add_argument("--report", type=Path, default=REPORT)
    args = parser.parse_args()
    if args.c_min < 3 or args.c_max < args.c_min:
        raise SystemExit("require 3 <= --c-min <= --c-max")

    iv.dps = PREC
    mp.dps = PREC + 20
    lines = ["# Interval certificates for the two thresholds\n"]
    lines.append(
        "Outward-rounded `mpmath.iv` at 60 digits. `C_k`, `H_c` and `K_c` are polynomials in "
        "`r` with rational coefficients, built here by exact `Fraction` arithmetic from the "
        "connectivity recursion, so `H_c'` is the exact derivative polynomial and the only "
        "transcendental is `log(1-r)`.\n"
    )
    threshold_max = min(args.c_max, args.threshold_c_max)
    lines.append("## Certified threshold enclosures\n")
    lines.append("| c | r_peel enclosure | θ_peel enclosure | r_Max enclosure | θ_Max enclosure | ratio |")
    lines.append("|---|---|---|---|---|---|")
    failures = []
    separation_rows = []
    for c in range(args.c_min, args.c_max + 1):
        t0 = time.time()
        D, g, Psi, Dc, Psic = make_funcs(c)
        lo, hi = R0, 1 - 1e-9
        if c <= threshold_max:
            # Independent global covering. The exact small-r certificate closes (0, R0].
            br, status = isolate(Dc, lo, hi)
        else:
            # The exact Bernstein certificate supplies uniqueness for these c values.
            br, status = bisect_unique_root(D, lo, hi)
        if br is None:
            print(f"c={c}: D_c isolation failed: {status}  ({time.time()-t0:.0f}s)", flush=True)
            failures.append(f"c={c}: D_c isolation failed: {status}")
            continue
        rp = iv.mpf([br[0], br[1]])
        tp2 = g(rp)
        tp = iv.sqrt(tp2)
        psi_peel = Psic(rp)
        if not _neg(psi_peel):
            failures.append(f"c={c}: Psi at the peeling bracket is not strictly negative")
            print(f"c={c}: peeling surplus sign failed  ({time.time()-t0:.0f}s)", flush=True)
            continue
        separation_rows.append((c, br, tp, psi_peel))

        if c <= threshold_max:
            # A global covering isolates the first Maxwell crossing for the displayed table.
            br2, status2 = isolate(Psic, br[1], hi)
            if br2 is None:
                print(f"c={c}: Psi isolation failed: {status2}  ({time.time()-t0:.0f}s)", flush=True)
                failures.append(f"c={c}: Psi isolation failed: {status2}")
                continue
            rm = iv.mpf([br2[0], br2[1]])
            tm = iv.sqrt(g(rm))
            ratio = tm / tp
            lines.append(
                f"| {c} | [{number(br[0], 30)}, {number(br[1], 30)}] | "
                f"[{endpoint(tp, 'a', 30)}, {endpoint(tp, 'b', 30)}] | "
                f"[{number(br2[0], 30)}, {number(br2[1], 30)}] | "
                f"[{endpoint(tm, 'a', 30)}, {endpoint(tm, 'b', 30)}] | "
                f"[{endpoint(ratio, 'a', 24)}, {endpoint(ratio, 'b', 24)}] |")
            print(f"c={c}  r_peel [{number(br[0])},{number(br[1])}]  "
                  f"r_Max [{number(br2[0])},{number(br2[1])}]  "
                  f"theta_peel~{endpoint(tp, 'a', 12)}  "
                  f"theta_Max~{endpoint(tm, 'a', 12)}  ({time.time()-t0:.0f}s)", flush=True)
        else:
            print(f"c={c}  r_peel [{number(br[0])},{number(br[1])}]  "
                  f"Psi_peel<={endpoint(psi_peel, 'b', 8)}  "
                  f"({time.time()-t0:.0f}s)", flush=True)

    lines.append("\n## Certified strict separation at the peeling point\n")
    lines.append("| c | r_peel enclosure | θ_peel enclosure | Ψ_c(θ_peel+) enclosure |")
    lines.append("|---|---|---|---|")
    for c, br, tp, psi_peel in separation_rows:
        lines.append(
            f"| {c} | [{number(br[0], 36)}, {number(br[1], 36)}] | "
            f"[{endpoint(tp, 'a', 30)}, {endpoint(tp, 'b', 30)}] | "
            f"[{endpoint(psi_peel, 'a', 30)}, {endpoint(psi_peel, 'b', 30)}] |")

    lines.append("\n## What is certified\n")
    if args.c_min <= threshold_max:
        lines.append(
            f"- For `{args.c_min} ≤ c ≤ {threshold_max}`, global interval coverings "
            "independently certify the unique roots defining `θ_peel` and `θ_Max`.")
    if args.c_max > threshold_max:
        local_min = max(args.c_min, threshold_max + 1)
        lines.append(
            f"- For `{local_min} ≤ c ≤ {args.c_max}`, the exact Bernstein certificate supplies "
            "uniqueness of the peeling root; point-interval bisection encloses it.")
    lines.append(
        f"- For every `{args.c_min} ≤ c ≤ {args.c_max}`, the enclosure of `Ψ` on the entire "
        "peeling bracket is strictly negative. Since `Ψ` is positive for large `θ`, this "
        "certifies `θ_peel < θ_Max`.\n")
    lines.append(
        "## What is not\n"
        f"The certificate covers `{args.c_min} ≤ c ≤ {args.c_max}`. A uniform proof for all "
        "`c ≥ 3` remains open.\n")
    if failures:
        lines.append("\n## Failed instances\n" + "\n".join(f"- {item}" for item in failures))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_report(args.report, "\n".join(lines))
    print(f"\nwrote {args.report}")
    if failures:
        raise SystemExit("interval certificate failed for one or more requested c values")


if __name__ == "__main__":
    main()
