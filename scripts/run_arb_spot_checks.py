"""Independent Arb spot checks for selected certificate instances.

This script does not import the mpmath interval implementation. It reconstructs the
connectivity recursion directly with python-flint/Arb ball arithmetic, reads only the decimal
root brackets emitted by ``run_interval_certificate.py``, and checks their endpoint signs and
the surplus sign on the whole peeling bracket. Decimal intervals in the report are obtained by
directed rounding of Arb's exact dyadic endpoints.
"""

from __future__ import annotations

import argparse
import re
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

from flint import arb, ctx, fmpq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from prfield.report import write_report  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
SOURCE_REPORT = ROOT / "reports" / "AUDIT_INTERVAL_CERTIFICATE.md"
REPORT = ROOT / "reports" / "AUDIT_ARB_SPOT_CHECK.md"
SPOT_CASES = (3, 4, 7, 24)
DECIMAL_PLACES = 28


class Dual:
    __slots__ = ("v", "d")

    def __init__(self, value: arb, derivative: arb):
        self.v = value
        self.d = derivative

    def __add__(self, other: object) -> "Dual":
        other = lift(other)
        return Dual(self.v + other.v, self.d + other.d)

    __radd__ = __add__

    def __sub__(self, other: object) -> "Dual":
        other = lift(other)
        return Dual(self.v - other.v, self.d - other.d)

    def __rsub__(self, other: object) -> "Dual":
        other = lift(other)
        return Dual(other.v - self.v, other.d - self.d)

    def __mul__(self, other: object) -> "Dual":
        other = lift(other)
        return Dual(self.v * other.v, self.d * other.v + self.v * other.d)

    __rmul__ = __mul__


def lift(value: object) -> Dual:
    return value if isinstance(value, Dual) else Dual(arb(value), arb(0))


def q_power(r: arb, exponent: int) -> Dual:
    q = arb(1) - r
    if exponent == 0:
        return Dual(arb(1), arb(0))
    return Dual(q**exponent, -exponent * q ** (exponent - 1))


def connected_duals(c: int, r: arb) -> dict[int, Dual]:
    connected = {1: Dual(arb(1), arb(0))}
    for size in range(2, c + 1):
        accumulator = Dual(arb(0), arb(0))
        for first_size in range(1, size):
            accumulator += (
                connected[first_size]
                * q_power(r, first_size * (size - first_size))
                * comb(size - 1, first_size - 1)
            )
        connected[size] = Dual(arb(1), arb(0)) - accumulator
    return connected


def evaluate(c: int, r: arb) -> tuple[arb, arb, arb]:
    connected = connected_duals(c, r)
    h = Dual(arb(0), arb(0))
    k_value = arb(0)
    for size in range(1, c + 1):
        factor = connected[size] * q_power(r, size * (c - size))
        h += factor * fmpq(comb(c - 1, size - 1) * (size - 1), c - 1)
        k_value += factor.v * comb(c, size)
    s = -(arb(1) - r).log()
    d_value = h.v - s * (arb(1) - r) * h.d
    psi = 2 * (k_value - c) + comb(c, 2) * s * (2 - h.v)
    theta = (s / h.v).sqrt()
    return d_value, psi, theta


def decimal_fraction(text: str) -> fmpq:
    value = Fraction(text)
    return fmpq(value.numerator, value.denominator)


def point(text: str) -> arb:
    return arb(decimal_fraction(text))


def interval(lo: str, hi: str) -> arb:
    lower, upper = decimal_fraction(lo), decimal_fraction(hi)
    if lower > upper:
        raise ValueError(f"reversed interval [{lo}, {hi}]")
    return arb((lower + upper) / 2, (upper - lower) / 2)


def parse_pair(cell: str) -> tuple[str, str]:
    match = re.fullmatch(r"\[\s*([^,]+),\s*([^]]+)\s*\]", cell.strip())
    if not match:
        raise ValueError(f"cannot parse interval cell: {cell}")
    return match.group(1).strip(), match.group(2).strip()


def read_brackets(path: Path) -> tuple[dict[int, tuple[str, str]], dict[int, tuple[str, str]]]:
    peeling: dict[int, tuple[str, str]] = {}
    maxwell: dict[int, tuple[str, str]] = {}
    section = ""
    for line in path.read_text().splitlines():
        if line.startswith("## "):
            section = line[3:]
            continue
        if not line.startswith("|") or line.startswith("|---"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if not cells or not cells[0].isdigit():
            continue
        c = int(cells[0])
        if section == "Certified threshold enclosures":
            peeling[c] = parse_pair(cells[1])
            maxwell[c] = parse_pair(cells[3])
        elif section == "Certified strict separation at the peeling point":
            peeling[c] = parse_pair(cells[1])
    return peeling, maxwell


def rounded_decimal(value: fmpq, places: int, upper: bool) -> str:
    numerator, denominator = int(value.p), int(value.q)
    scale = 10**places
    scaled = -((-numerator * scale) // denominator) if upper else (numerator * scale) // denominator
    sign = "-" if scaled < 0 else ""
    digits = str(abs(scaled)).rjust(places + 1, "0")
    return f"{sign}{digits[:-places]}.{digits[-places:]}"


def enclosure(value: arb, places: int = DECIMAL_PLACES) -> str:
    return (
        f"[{rounded_decimal(value.lower().fmpq(), places, False)}, "
        f"{rounded_decimal(value.upper().fmpq(), places, True)}]"
    )


def strictly_negative(value: arb) -> bool:
    return value.upper() < 0


def strictly_positive(value: arb) -> bool:
    return value.lower() > 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-report", type=Path, default=SOURCE_REPORT)
    parser.add_argument("--report", type=Path, default=REPORT)
    args = parser.parse_args()

    ctx.dps = 100
    peeling, maxwell = read_brackets(args.source_report)
    failures: list[str] = []
    rows: list[str] = []
    for c in SPOT_CASES:
        if c not in peeling:
            failures.append(f"c={c}: peeling bracket missing")
            continue
        peel_lo, peel_hi = peeling[c]
        d_lo, _, _ = evaluate(c, point(peel_lo))
        d_hi, _, _ = evaluate(c, point(peel_hi))
        _, psi_peel, theta_peel = evaluate(c, interval(peel_lo, peel_hi))
        peel_ok = strictly_negative(d_lo) and strictly_positive(d_hi) and strictly_negative(psi_peel)

        maxwell_result = "not tabulated"
        if c in maxwell:
            max_lo, max_hi = maxwell[c]
            _, psi_lo, _ = evaluate(c, point(max_lo))
            _, psi_hi, _ = evaluate(c, point(max_hi))
            _, psi_max, theta_max = evaluate(c, interval(max_lo, max_hi))
            max_ok = strictly_negative(psi_lo) and strictly_positive(psi_hi) and 0 in psi_max
            maxwell_result = (
                f"{enclosure(theta_max)}; Psi(root box)={enclosure(psi_max)}; "
                f"{'PASS' if max_ok else 'FAIL'}"
            )
            if not max_ok:
                failures.append(f"c={c}: Arb does not verify the Maxwell bracket")

        rows.append(
            f"| {c} | {enclosure(theta_peel)} | {enclosure(psi_peel)} | "
            f"{'PASS' if peel_ok else 'FAIL'} | {maxwell_result} |"
        )
        if not peel_ok:
            failures.append(f"c={c}: Arb does not verify the peeling bracket and surplus sign")

    lines = [
        "# Independent Arb spot checks\n",
        "Generated by `scripts/run_arb_spot_checks.py` with `python-flint`/Arb at 100 decimal "
        "digits. The implementation reconstructs the connectivity recursion independently of "
        "the `mpmath.iv` code. Every displayed endpoint is obtained by directed decimal "
        "rounding of an exact dyadic Arb endpoint.\n",
        "| c | Arb theta_peel enclosure | Arb Psi on peeling box | peeling check | Maxwell check |",
        "|---|---|---|---|---|",
        *rows,
        "\nThe peeling check requires `D(lo)<0<D(hi)` and `Psi<0` on the complete peeling "
        "box. Where a Maxwell bracket is tabulated, the check requires opposite endpoint signs "
        "and an Arb enclosure of `Psi` on the root box that contains zero.\n",
    ]
    if failures:
        lines += ["## Failures", *(f"- {failure}" for failure in failures)]
    write_report(args.report, "\n".join(lines))
    print(f"wrote {args.report.relative_to(ROOT)}")
    for row in rows:
        print(row)
    if failures:
        raise SystemExit("Arb spot check failed")


if __name__ == "__main__":
    main()
