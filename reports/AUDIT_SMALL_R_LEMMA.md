<!-- GENERATED: everything above the ANALYSIS marker is written by the script named in
     this file's first paragraph and is overwritten on every rerun. -->

# The small-`r` lemma

`D_c(0) = 0` to second order, so no interval covering can certify a sign near the origin. This closes `(0, r0]` by an exact series argument instead, with `r0 = 0.005` and truncation order `N = 48`.

| c | d₀ | d₁ | d₂ (exact) | −(2c−5)/2 | ‖H_c'‖₁ | bound on D_c/r² | certified |
|---|---|---|---|---|---|---|---|
| 3 | 0 | 0 | -1/2 | -1/2 | 6 | -0.484140 | yes |
| 4 | 0 | 0 | -3/2 | -3/2 | 80 | -1.488618 | yes |
| 5 | 0 | 0 | -5/2 | -5/2 | 3.6e+03 | -2.484558 | yes |
| 6 | 0 | 0 | -7/2 | -7/2 | 4.18e+05 | -3.439088 | yes |
| 7 | 0 | 0 | -9/2 | -9/2 | 1.19e+08 | -4.374333 | yes |
| 8 | 0 | 0 | -11/2 | -11/2 | 7.91e+10 | -5.290402 | yes |
| 9 | 0 | 0 | -13/2 | -13/2 | 1.19e+14 | -6.182100 | yes |
| 10 | 0 | 0 | -15/2 | -15/2 | 4.02e+17 | -7.051724 | yes |
| 11 | 0 | 0 | -17/2 | -17/2 | 2.99e+21 | -7.898868 | yes |
| 12 | 0 | 0 | -19/2 | -19/2 | 4.88e+25 | -8.723139 | yes |
| 13 | 0 | 0 | -21/2 | -21/2 | 1.73e+30 | -9.524091 | yes |
| 14 | 0 | 0 | -23/2 | -23/2 | 1.32e+35 | -10.300230 | yes |
| 15 | 0 | 0 | -25/2 | -25/2 | 2.16e+40 | -11.051635 | yes |
| 16 | 0 | 0 | -27/2 | -27/2 | 7.56e+45 | -11.777756 | yes |
| 17 | 0 | 0 | -29/2 | -29/2 | 5.63e+51 | -12.478037 | yes |
| 18 | 0 | 0 | -31/2 | -31/2 | 8.88e+57 | -13.151910 | yes |
| 19 | 0 | 0 | -33/2 | -33/2 | 2.96e+64 | -13.798801 | yes |
| 20 | 0 | 0 | -35/2 | -35/2 | 2.08e+71 | -14.418026 | yes |
| 21 | 0 | 0 | -37/2 | -37/2 | 3.06e+78 | -15.008656 | yes |
| 22 | 0 | 0 | -39/2 | -39/2 | 9.46e+85 | -15.570264 | yes |
| 23 | 0 | 0 | -41/2 | -41/2 | 6.11e+93 | -16.102177 | yes |
| 24 | 0 | 0 | -43/2 | -43/2 | 8.25e+101 | -16.603702 | yes |

## What this establishes
`d₀ = d₁ = 0` and `d₂ = −(2c−5)/2` exactly, so the vanishing at the origin is exactly second order and the leading coefficient is at most `−1/2` for every `c ≥ 3`. The bound column is a single exact rational upper bound for `D_c(r)/r²` on the whole of `(0, r0]`, computed from the exact Taylor coefficients plus the explicit tail bound; being negative, it certifies `D_c < 0` there.

For the globally covered cases this closes the interval certificate at the origin. Across the full tested range it independently confirms the left-endpoint sign used by the exact Bernstein single-crossing certificate.

## What it does not
This is a per-instance certificate for `3 ≤ c ≤ 24`. A uniform proof for all `c ≥ 3` requires a bound on the Taylor coefficients that is independent of `c`.


<!-- ANALYSIS: hand-written below; preserved across reruns -->
