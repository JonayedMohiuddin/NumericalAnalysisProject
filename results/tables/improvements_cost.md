**Mean LU factorizations on the cases solved by both this combination and the paper's method (ours vs paper), and median run time relative to the paper's method**

| options | S1 paper | S2 Sec 4.4 | S3 Sec 4.4 | S4 weak K | S5 strong K | S1 median time ratio |
|---|---|---|---|---|---|---|
| paper | 8.30 vs 8.30 | 8.59 vs 8.59 | 8.68 vs 8.68 | 8.59 vs 8.59 | 8.22 vs 8.22 | 1.00 |
| adaptive | 8.83 vs 8.30 | 8.79 vs 8.59 | 8.75 vs 8.68 | 9.00 vs 8.59 | 8.52 vs 8.13 | 0.99 |
| PC | 8.60 vs 8.30 | 8.69 vs 8.59 | 8.64 vs 8.68 | 9.23 vs 8.59 | 9.00 vs 8.22 | 1.00 |
| PC + adaptive | 8.87 vs 8.30 | 8.79 vs 8.59 | 8.64 vs 8.68 | 10.62 vs 8.71 | 19.03 vs 8.22 | 1.01 |
| OM | 7.37 vs 8.30 | 7.55 vs 8.59 | 7.61 vs 8.68 | 7.55 vs 8.59 | 7.34 vs 8.22 | 0.92 |
| OM + adaptive | 7.90 vs 8.30 | 7.76 vs 8.59 | 7.68 vs 8.68 | 7.95 vs 8.59 | 7.68 vs 8.13 | 0.92 |
| OM + PC | 7.70 vs 8.30 | 7.79 vs 8.59 | 7.71 vs 8.68 | 8.05 vs 8.59 | 8.06 vs 8.22 | 0.98 |
| OM + PC + adaptive | 7.97 vs 8.30 | 7.90 vs 8.59 | 7.71 vs 8.68 | 9.43 vs 8.71 | 18.50 vs 8.22 | 0.98 |
| BE-chord | 8.75 vs 7.88 | 8.83 vs 7.96 | 8.67 vs 7.86 | 8.67 vs 7.80 | 8.64 vs 7.89 | 1.06 |
| RK4 | 7.00 vs 3.00 | 7.00 vs 3.33 | 9.60 vs 5.20 | 7.00 vs 3.00 | 10.00 vs 3.00 | 2.30 |
| scaled homotopy | 8.30 vs 8.30 | 8.26 vs 8.33 | 8.54 vs 8.50 | 8.59 vs 8.59 | 8.31 vs 8.22 | 1.04 |
| newton homotopy | 8.94 vs 8.00 | 8.83 vs 8.11 | 9.38 vs 8.12 | 8.94 vs 7.88 | 8.94 vs 8.35 | 1.10 |
