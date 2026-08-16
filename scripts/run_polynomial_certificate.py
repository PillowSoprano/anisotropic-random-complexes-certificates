"""Independent exact-polynomial check of the shape condition.

The interval certificate checks D_c directly.  This script follows the independent route
stated in Appendix P2: with W(t)=H_c(1-exp(-t)),

    sign W''(t) = sign G_c(r),  G_c(r)=(1-r)H_c''(r)-H_c'(r).

It reconstructs H_c from the connectivity recursion without importing any certificate code,
then uses a rational Sturm sequence to count the roots of G_c in (0,1).  A single root, with
G_c positive near 0 and negative near 1, makes W convex then concave.  Therefore
F(t)=W(t)-tW'(t) decreases from 0 and then increases to 1, so F (and hence D_c) changes sign
exactly once, negative then positive.

Usage:  PYTHONPATH=src python3 scripts/run_polynomial_certificate.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "AUDIT_POLYNOMIAL_CERTIFICATE.md"


def trim(p: list[F]) -> list[F]:
    p = p[:]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [F(0)]


def add(a: list[F], b: list[F]) -> list[F]:
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else F(0)) +
                 (b[i] if i < len(b) else F(0)) for i in range(n)])


def scale(a: list[F], s: F) -> list[F]:
    return trim([s * x for x in a])


def mul(a: list[F], b: list[F]) -> list[F]:
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def deriv(a: list[F]) -> list[F]:
    return trim([F(i) * a[i] for i in range(1, len(a))] or [F(0)])


def evaluate(a: list[F], x: F) -> F:
    out = F(0)
    for coefficient in reversed(a):
        out = out * x + coefficient
    return out


def divmod_poly(a: list[F], b: list[F]) -> tuple[list[F], list[F]]:
    a, b = trim(a), trim(b)
    if b == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(a) < len(b):
        return [F(0)], a
    q = [F(0)] * (len(a) - len(b) + 1)
    r = a[:]
    while r != [0] and len(r) >= len(b):
        shift = len(r) - len(b)
        coefficient = r[-1] / b[-1]
        q[shift] += coefficient
        for j, value in enumerate(b):
            r[j + shift] -= coefficient * value
        r = trim(r)
    return trim(q), trim(r)


def monic(a: list[F]) -> list[F]:
    a = trim(a)
    return scale(a, F(1, 1) / a[-1]) if a != [0] else a


def gcd_poly(a: list[F], b: list[F]) -> list[F]:
    a, b = trim(a), trim(b)
    while b != [0]:
        _, r = divmod_poly(a, b)
        a, b = b, r
    return monic(a)


def one_minus_r_pow(k: int) -> list[F]:
    return [F((-1) ** i * comb(k, i)) for i in range(k + 1)]


def connected_polys(c: int) -> dict[int, list[F]]:
    connected = {1: [F(1)]}
    for m in range(2, c + 1):
        acc = [F(0)]
        for k in range(1, m):
            term = mul(connected[k], one_minus_r_pow(k * (m - k)))
            acc = add(acc, scale(term, F(comb(m - 1, k - 1))))
        connected[m] = add([F(1)], scale(acc, F(-1)))
    return connected


def h_poly(c: int) -> list[F]:
    connected = connected_polys(c)
    out = [F(0)]
    for k in range(1, c + 1):
        weight = F(comb(c - 1, k - 1) * (k - 1), c - 1)
        out = add(out, scale(mul(connected[k], one_minus_r_pow(k * (c - k))), weight))
    return trim(out)


def g_poly(c: int) -> list[F]:
    hp = deriv(h_poly(c))
    hpp = deriv(hp)
    return add(mul([F(1), F(-1)], hpp), scale(hp, F(-1)))


def remove_one_endpoint_factors(p: list[F]) -> tuple[list[F], int]:
    multiplicity = 0
    factor = [F(1), F(-1)]
    while evaluate(p, F(1)) == 0:
        q, r = divmod_poly(p, factor)
        if r != [0]:
            raise AssertionError("failed to remove exact (1-r) factor")
        p = q
        multiplicity += 1
    return p, multiplicity


def sturm_sequence(p: list[F]) -> list[list[F]]:
    common = gcd_poly(p, deriv(p))
    square_free, remainder = divmod_poly(p, common)
    if remainder != [0]:
        raise AssertionError("square-free division was not exact")
    sequence = [monic(square_free), deriv(monic(square_free))]
    while sequence[-1] != [0]:
        _, r = divmod_poly(sequence[-2], sequence[-1])
        if r == [0]:
            break
        sequence.append(scale(r, F(-1)))
    return sequence


def variations(sequence: list[list[F]], x: F) -> int:
    signs: list[int] = []
    for p in sequence:
        value = evaluate(p, x)
        if value:
            signs.append(1 if value > 0 else -1)
    return sum(a != b for a, b in zip(signs, signs[1:]))


def certify(c: int) -> dict[str, int | bool]:
    g = g_poly(c)
    reduced, endpoint_multiplicity = remove_one_endpoint_factors(g)
    sequence = sturm_sequence(reduced)
    roots = variations(sequence, F(0)) - variations(sequence, F(1))
    g0 = evaluate(g, F(0))
    near_one_sign = evaluate(reduced, F(1))
    certified = roots == 1 and g0 == 2 * c - 5 and near_one_sign < 0
    return {
        "c": c,
        "degree": len(g) - 1,
        "endpoint_multiplicity": endpoint_multiplicity,
        "roots": roots,
        "g0": int(g0),
        "near_one_negative": near_one_sign < 0,
        "certified": certified,
    }


def main() -> None:
    # Cheap reconstruction check, independent of the interval implementation.
    assert h_poly(3) == [F(0), F(1), F(1), F(-1)]
    rows = [certify(c) for c in range(3, 10)]
    if not all(bool(row["certified"]) for row in rows):
        raise SystemExit("polynomial certificate failed")

    lines = [
        "# Independent exact-polynomial certificate\n",
        "Generated by `scripts/run_polynomial_certificate.py`. The script reconstructs "
        "`H_c` independently, removes endpoint factors `(1-r)` exactly, and applies Sturm's "
        "theorem over `Fraction` coefficients. No floating-point or interval evaluation is "
        "used.\n",
        "| c | deg G_c | multiplicity at r=1 | roots in (0,1) | G_c(0) | sign near 1 | status |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['c']} | {row['degree']} | {row['endpoint_multiplicity']} | "
            f"{row['roots']} | {row['g0']} | negative | CERTIFIED |"
        )
    lines += [
        "\nFor every `3 <= c <= 9`, `G_c` has exactly one interior root, is positive near "
        "zero and negative near one. Hence `W` is convex then concave, so "
        "`F(t)=W(t)-tW'(t)` and `D_c` change sign exactly once, negative then positive. This "
        "independently verifies the shape condition checked by the interval certificate.\n"
    ]
    REPORT.write_text("\n".join(lines))
    for row in rows:
        print(f"c={row['c']}  roots={row['roots']}  deg={row['degree']}  CERTIFIED")
    print(f"wrote {REPORT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
