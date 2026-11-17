## Per-scenario results (50 trials per cell)

| Scenario | Agents | M-MTS success | FE-CBF success | STCS success | M-MTS margin (m) | FE-CBF margin (m) | STCS margin (m) |
|---|---|---|---|---|---|---|---|
| F1 Opposing lanes | 3 | 50/50 | 0/50 | 0/50 | 1.69 ± 0.00 | 6.00 ± 0.00 | 6.00 ± 0.00 |
| F2 Perpendicular streams | 4 | 50/50 | 50/50 | 50/50 | 5.03 ± 0.16 | 3.04 ± 0.04 | 1.98 ± 0.19 |
| F3 Circle swap | 4 | 50/50 | 0/50 | 0/50 | 1.76 ± 0.02 | 21.47 ± 0.42 | 21.47 ± 0.42 |
| S1 Narrow Gate | 3 | 44/50 | 0/50 | 50/50 | 0.67 ± 0.03 | 0.75 ± 0.04 | 0.67 ± 0.03 |
| S2 Obstacle field | 4 | 48/50 | 0/50 | 50/50 | 0.73 ± 0.03 | 0.87 ± 0.06 | 0.72 ± 0.02 |
| S3 Walled corridor | 3 | 50/50 | 0/50 | 0/50 | 0.90 ± 0.01 | 4.90 ± 0.23 | 4.90 ± 0.23 |
| D1 Crossing obstacles | 2 | 50/50 | 48/50 (2 coll.) | 50/50 | 2.85 ± 0.26 | 2.57 ± 0.32 | 1.73 ± 0.08 |
| D2 Junction with obstacle | 3 | 50/50 | 50/50 | 50/50 | 3.88 ± 0.36 | 2.68 ± 0.14 | 1.79 ± 0.11 |
| D3 Obstacle stream | 3 | 46/50 (4 coll.) | 47/50 (3 coll.) | 49/50 (1 coll.) | 2.25 ± 0.27 | 2.09 ± 0.25 | 1.53 ± 0.11 |
| N1 Aggressive Cut-off | 2 | 50/50 | 50/50 | 49/50 (1 coll.) | 3.60 ± 0.56 | 2.77 ± 0.21 | 1.37 ± 0.11 |
| N2 Red-light runner | 2 | 50/50 | 50/50 | 50/50 | 5.27 ± 0.23 | 3.45 ± 0.23 | 2.16 ± 0.28 |
| N3 Weaving vehicle | 2 | 50/50 | 50/50 | 50/50 | 2.86 ± 0.15 | 2.80 ± 0.16 | 2.01 ± 0.19 |

## Table I: S1 and N1

| Scenario | Metric | M-MTS | FE-CBF | STCS |
|---|---|---|---|---|
| S1 | Success | 44/50 | 0/50 | 50/50 |
| S1 | Runtime (ms) | 26.80 ± 4.56 | 18.69 ± 3.43 | 2.40 ± 0.16 |
| S1 | Safety margin (m) | 0.67 ± 0.03 | 0.75 ± 0.04 | 0.67 ± 0.03 |
| S1 | Suboptimality | 2.42 ± 0.07 | n/a | 1.27 ± 0.02 |
| S1 | Tracking error (m) | 0.078 ± 0.012 | n/a | 0.043 ± 0.002 |
| Scenario | Metric | M-MTS | FE-CBF | STCS |
|---|---|---|---|---|
| N1 | Success | 50/50 | 50/50 | 49/50 |
| N1 | Runtime (ms) | 15.28 ± 0.34 | 4.58 ± 0.07 | 5.21 ± 0.39 |
| N1 | Safety margin (m) | 3.60 ± 0.56 | 2.77 ± 0.21 | 1.37 ± 0.11 |
| N1 | Suboptimality | 2.19 ± 0.02 | 2.05 ± 0.01 | 1.94 ± 0.01 |
| N1 | Tracking error (m) | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.000 ± 0.000 |

## Statistical tests

- F1, M-MTS vs FE-CBF: success Fisher p = 1.98e-29; margin Welch t = -1884.63 (df ≈ 49), p = 9.51e-121; suboptimality t = nan, p = nan; runtime t = -151.21, p = 2.36e-68
- F1, M-MTS vs STCS: success Fisher p = 1.98e-29; margin Welch t = -1884.63 (df ≈ 49), p = 9.51e-121; suboptimality t = nan, p = nan; runtime t = -111.34, p = 1.34e-63
- F2, M-MTS vs FE-CBF: success Fisher p = 1; margin Welch t = 24.51 (df ≈ 57), p = 6.31e-32; suboptimality t = 9.57, p = 1.01e-13; runtime t = 30.56, p = 1.71e-34
- F2, M-MTS vs STCS: success Fisher p = 1; margin Welch t = 24.72 (df ≈ 94), p = 5.36e-43; suboptimality t = 30.78, p = 1.01e-35; runtime t = 37.34, p = 1.43e-38
- F3, M-MTS vs FE-CBF: success Fisher p = 1.98e-29; margin Welch t = -94.65 (df ≈ 49), p = 2.57e-57; suboptimality t = nan, p = nan; runtime t = -121.88, p = 3.22e-85
- F3, M-MTS vs STCS: success Fisher p = 1.98e-29; margin Welch t = -94.65 (df ≈ 49), p = 2.57e-57; suboptimality t = nan, p = nan; runtime t = -50.90, p = 1.87e-52
- S1, M-MTS vs FE-CBF: success Fisher p = 6.44e-22; margin Welch t = -3.06 (df ≈ 90), p = 0.00291; suboptimality t = nan, p = nan; runtime t = 2.86, p = 0.00526
- S1, M-MTS vs STCS: success Fisher p = 0.0267; margin Welch t = -0.16 (df ≈ 96), p = 0.87; suboptimality t = 31.53, p = 1.82e-33; runtime t = 10.75, p = 1.65e-14
- S2, M-MTS vs FE-CBF: success Fisher p = 2.63e-26; margin Welch t = -4.38 (df ≈ 69), p = 4.16e-05; suboptimality t = nan, p = nan; runtime t = 5.80, p = 1.03e-07
- S2, M-MTS vs STCS: success Fisher p = 0.495; margin Welch t = 0.67 (df ≈ 97), p = 0.504; suboptimality t = 46.06, p = 6.76e-43; runtime t = 14.57, p = 1.86e-19
- S3, M-MTS vs FE-CBF: success Fisher p = 1.98e-29; margin Welch t = -34.60 (df ≈ 49), p = 3.42e-36; suboptimality t = nan, p = nan; runtime t = -125.49, p = 5.36e-65
- S3, M-MTS vs STCS: success Fisher p = 1.98e-29; margin Welch t = -34.60 (df ≈ 49), p = 3.42e-36; suboptimality t = nan, p = nan; runtime t = -98.20, p = 8.16e-64
- D1, M-MTS vs FE-CBF: success Fisher p = 0.495; margin Welch t = 1.34 (df ≈ 94), p = 0.183; suboptimality t = 4.42, p = 2.6e-05; runtime t = 11.06, p = 2.54e-18
- D1, M-MTS vs STCS: success Fisher p = 1; margin Welch t = 8.19 (df ≈ 59), p = 2.72e-11; suboptimality t = 16.59, p = 1.6e-28; runtime t = 27.65, p = 3.17e-32
- D2, M-MTS vs FE-CBF: success Fisher p = 1; margin Welch t = 6.21 (df ≈ 63), p = 4.78e-08; suboptimality t = 11.65, p = 3.43e-20; runtime t = 30.01, p = 1.59e-36
- D2, M-MTS vs STCS: success Fisher p = 1; margin Welch t = 10.99 (df ≈ 59), p = 7.18e-16; suboptimality t = 35.03, p = 3.86e-56; runtime t = 44.57, p = 7.67e-43
- D3, M-MTS vs FE-CBF: success Fisher p = 1; margin Welch t = 0.87 (df ≈ 97), p = 0.386; suboptimality t = 1.27, p = 0.207; runtime t = 5.07, p = 5.71e-06
- D3, M-MTS vs STCS: success Fisher p = 0.362; margin Welch t = 4.94 (df ≈ 63), p = 6e-06; suboptimality t = 5.30, p = 1.38e-06; runtime t = 6.06, p = 1.87e-07
- N1, M-MTS vs FE-CBF: success Fisher p = 1; margin Welch t = 2.80 (df ≈ 62), p = 0.00687; suboptimality t = 13.19, p = 4.47e-23; runtime t = 62.19, p = 4.48e-51
- N1, M-MTS vs STCS: success Fisher p = 1; margin Welch t = 7.86 (df ≈ 53), p = 1.86e-10; suboptimality t = 24.60, p = 1.06e-39; runtime t = 39.23, p = 5.12e-61
- N2, M-MTS vs FE-CBF: success Fisher p = 1; margin Welch t = 11.16 (df ≈ 98), p = 3.78e-19; suboptimality t = 6.33, p = 7.64e-09; runtime t = 31.67, p = 5.85e-41
- N2, M-MTS vs STCS: success Fisher p = 1; margin Welch t = 17.13 (df ≈ 95), p = 8.61e-31; suboptimality t = 35.39, p = 1.29e-55; runtime t = 45.42, p = 5.41e-51
- N3, M-MTS vs FE-CBF: success Fisher p = 1; margin Welch t = 0.50 (df ≈ 98), p = 0.619; suboptimality t = 6.33, p = 7.56e-09; runtime t = 75.18, p = 3.86e-64
- N3, M-MTS vs STCS: success Fisher p = 1; margin Welch t = 7.01 (df ≈ 93), p = 3.69e-10; suboptimality t = 56.68, p = 1.77e-75; runtime t = 89.62, p = 1.13e-66

## Minimum clearance over trials in which the agents drove (success or collision; time-outs excluded)

| Percentile | M-MTS | FE-CBF | STCS |
|---|---|---|---|
| 5 | 0.66 | 1.56 | 0.61 |
| 50 | 2.25 | 2.83 | 1.54 |
| 95 | 5.73 | 4.09 | 2.87 |
| min | -0.14 | -0.22 | -0.02 |
| collisions | 4/592 | 5/350 | 2/450 |

## Overall

- M-MTS: {"success": 588, "n": 600, "collisions": 4, "timeouts": 8, "qp_inf": 0.004239233352268211, "dfcs_fail": 0.0023808080808080806, "dfcs_budget": 0.0, "track": 0.03909164844272334, "track_max": 0.40070199977530213, "runtime": 21.56387023305724, "runtime_max": 529.3408000006821, "qp_ms": 8.625979104129032}
- FE-CBF: {"success": 345, "n": 600, "collisions": 5, "timeouts": 250, "qp_inf": 0.08655587212467962, "dfcs_fail": 0.25606881313131313, "dfcs_budget": 0.25243750000000004, "track": 0.003988776693361008, "track_max": 0.09593564024477086, "runtime": 39.016047467984116, "runtime_max": 203.47840000067663, "qp_ms": 0.334248149596168}
- STCS: {"success": 448, "n": 600, "collisions": 2, "timeouts": 150, "qp_inf": 0.0, "dfcs_fail": 0.2512003532182103, "dfcs_budget": 0.25036309523809525, "track": 0.009315171909896976, "track_max": 0.056506410080650486, "runtime": 18.85388189921851, "runtime_max": 154.32789999977103, "qp_ms": NaN}

## Scaling world

| N | M-MTS success | FE-CBF success | STCS success | M-MTS runtime (ms) | FE-CBF runtime (ms) | STCS runtime (ms) |
|---|---|---|---|---|---|---|
| 2 | 30/30 | 6/30 | 30/30 | 6.0 ± 0.7 | 4.7 ± 1.6 | 1.5 ± 0.2 |
| 4 | 30/30 | 3/30 | 30/30 | 17.7 ± 2.0 | 18.2 ± 9.5 | 2.7 ± 0.3 |
| 6 | 30/30 | 1/30 | 30/30 | 45.7 ± 6.2 | 26.2 ± 7.8 | 4.1 ± 0.3 |
| 8 | 30/30 | 0/30 | 30/30 | 73.0 ± 6.7 | 40.1 ± 6.7 | 5.7 ± 0.4 |
| 10 | 30/30 | 0/30 | 30/30 | 125.0 ± 15.6 | 52.9 ± 8.3 | 8.2 ± 0.7 |
| 12 | 30/30 | 0/30 | 30/30 | 212.7 ± 27.2 | 80.4 ± 10.7 | 9.0 ± 0.7 |
| 14 | 30/30 | 0/30 | 30/30 | 359.7 ± 43.2 | 105.8 ± 14.1 | 13.0 ± 0.9 |
| 16 | 29/30 | 0/30 | 30/30 | 473.2 ± 57.2 | 134.0 ± 15.6 | 15.9 ± 1.1 |
