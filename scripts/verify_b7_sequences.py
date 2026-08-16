#!/usr/bin/env python3
"""Exact checks for the Bernstein reduction in Appendix B.7.

The script has no third-party dependencies.  It constructs the two-terminal
reliability polynomial H_c from the exact connected-graph recursion, converts
H_c to degree-M Bernstein form, and derives

    q_m = (M-m)(a_{m+1}-a_m),
    eta_m = (M-m)(q_{m+1}-q_m)/q_m.

All decisions use fractions.  Decimal values appear only in the descriptive
least-squares fit for the width of the central log-concave window.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb

F = Fraction


def trim(a: list[F]) -> list[F]:
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a: list[F], b: list[F]) -> list[F]:
    out = [F(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return trim(out)


def scale(a: list[F], value: F) -> list[F]:
    return trim([value * x for x in a])


def multiply(a: list[F], b: list[F]) -> list[F]:
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def one_minus_r_power(k: int) -> list[F]:
    return [F((-1) ** i * comb(k, i)) for i in range(k + 1)]


def connected_polynomials(c: int) -> dict[int, list[F]]:
    connected = {1: [F(1)]}
    for size in range(2, c + 1):
        accumulator = [F(0)]
        for first_size in range(1, size):
            term = multiply(
                connected[first_size],
                one_minus_r_power(first_size * (size - first_size)),
            )
            accumulator = add(
                accumulator,
                scale(term, F(comb(size - 1, first_size - 1))),
            )
        connected[size] = add([F(1)], scale(accumulator, F(-1)))
    return connected


def h_polynomial(c: int) -> list[F]:
    connected = connected_polynomials(c)
    out = [F(0)]
    for size in range(1, c + 1):
        weight = F(comb(c - 1, size - 1) * (size - 1), c - 1)
        term = multiply(
            connected[size], one_minus_r_power(size * (c - size))
        )
        out = add(out, scale(term, weight))
    return trim(out)


def power_to_bernstein(power: list[F], degree: int) -> list[F]:
    padded = power + [F(0)] * (degree + 1 - len(power))
    return [
        sum(
            padded[k] * F(comb(m, k), comb(degree, k))
            for k in range(m + 1)
        )
        for m in range(degree + 1)
    ]


def b7_sequences(c: int) -> tuple[list[F], list[F], list[F]]:
    edge_count = comb(c, 2)
    a = power_to_bernstein(h_polynomial(c), edge_count)
    q = [
        F(edge_count - m) * (a[m + 1] - a[m])
        for m in range(edge_count)
    ] + [F(0)]
    support_end = comb(c - 1, 2)
    eta = [
        F(edge_count - m) * (q[m + 1] - q[m]) / q[m]
        for m in range(support_end + 1)
    ]
    return a, q, eta


def sign_changes(values: list[F]) -> int:
    signs = [1 if value > 0 else -1 for value in values if value != 0]
    return sum(signs[i] != signs[i - 1] for i in range(1, len(signs)))


def central_window(q: list[F], support_end: int, mode: int) -> tuple[int, int]:
    if support_end < 2:
        return 0, 0
    failures = [
        m
        for m in range(1, support_end)
        if q[m] * q[m] < q[m - 1] * q[m + 1]
    ]
    left_failure = max((m for m in failures if m < mode), default=0)
    right_failure = min(
        (m for m in failures if m > mode), default=support_end
    )
    left = left_failure + 1 if left_failure else 1
    right = right_failure - 1 if right_failure < support_end else support_end - 1
    return left, right


def linear_fit(xs: list[int], ys: list[int]) -> tuple[float, float, float]:
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    denominator = sum((x - x_mean) ** 2 for x in xs)
    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denominator
    intercept = y_mean - slope * x_mean
    residual = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    total = sum((y - y_mean) ** 2 for y in ys)
    return slope, intercept, 1.0 - residual / total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-max", type=int, default=24)
    args = parser.parse_args()
    if args.c_max < 3:
        raise SystemExit("--c-max must be at least 3")

    fit_c: list[int] = []
    fit_gap: list[int] = []
    print(" c  mode  eta decreasing  eta recross  dq changes  LC window  right gap")
    for c in range(3, args.c_max + 1):
        edge_count = comb(c, 2)
        support_end = comb(c - 1, 2)
        _, q, eta = b7_sequences(c)

        assert q[0] == 1
        assert q[1] - q[0] == F(2 * c - 5, edge_count)
        assert q[2] - q[1] == F(
            (c - 2) * (2 * c * c + 9 * c - 61),
            2 * edge_count * (edge_count - 1),
        )
        assert all(q[m] > 0 for m in range(support_end + 1))
        assert all(q[m] == 0 for m in range(support_end + 1, edge_count + 1))

        mode = max(range(support_end + 1), key=lambda m: q[m])
        differences = [q[m + 1] - q[m] for m in range(support_end + 1)]
        difference_changes = sign_changes(differences)
        first_negative = next(m for m, value in enumerate(eta) if value < 0)
        recrosses = sum(value >= 0 for value in eta[first_negative:])
        eta_decreasing = all(eta[m + 1] <= eta[m] for m in range(support_end))
        left, right = central_window(q, support_end, mode)

        assert difference_changes == 1
        assert recrosses == 0
        if c == 4:
            assert eta == [F(3), F(7, 9), F(-43, 13), F(-3)]
        if c == 8:
            assert q[0] == 1 and q[1] == F(39, 28) and q[2] == F(35, 18)
            assert q[1] * q[1] - q[0] * q[2] == F(-31, 7056)
        if c >= 8:
            assert q[1] * q[1] < q[0] * q[2]

        if c >= 8:
            fit_c.append(c)
            fit_gap.append(right - mode)
        window = "empty" if support_end < 2 else f"[{left:2d},{right:3d}]"
        print(
            f"{c:2d}  {mode:4d}  {str(eta_decreasing):14s}"
            f"  {recrosses:11d}  {difference_changes:10d}"
            f"  {window:8s}  {right - mode:9d}"
        )

    if len(fit_c) >= 2:
        slope, intercept, r_squared = linear_fit(fit_c, fit_gap)
        print()
        print(
            "Descriptive fit for right-window excess: "
            f"r_c-mode = {slope:.6f} c + {intercept:.6f}, R^2={r_squared:.6f}"
        )


if __name__ == "__main__":
    main()
