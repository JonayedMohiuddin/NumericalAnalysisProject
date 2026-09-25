**Runs that reach the reference operating point, out of 34 per setting. 'other root' counts runs that converged to a different solution**

| options | S1 paper | S2 Sec 4.4 | S3 Sec 4.4 | S4 weak K | S5 strong K | total solved | other root | lost vs paper |
|---|---|---|---|---|---|---|---|---|
| paper | 28 | 27 | 25 | 18 | 29 | 127 | 14 | 0 |
| adaptive | 30 | 30 | 26 | 18 | 28 | 132 | 14 | 1 |
| PC | 32 | 32 | 32 | 22 | 34 | 152 | 2 | 0 |
| PC + adaptive | 32 | 32 | 32 | 22 | 34 | 152 | 1 | 0 |
| OM | 30 | 30 | 29 | 20 | 34 | 143 | 12 | 0 |
| OM + adaptive | 32 | 33 | 30 | 20 | 33 | 148 | 10 | 1 |
| OM + PC | 32 | 32 | 32 | 22 | 34 | 152 | 3 | 0 |
| OM + PC + adaptive | 32 | 32 | 32 | 22 | 34 | 152 | 2 | 0 |
| richardson | 32 | 34 | 28 | 23 | 34 | 151 | 9 | 1 |
| richardson + FE | 34 | 33 | 29 | 27 | 34 | 157 | 3 | 1 |
| richardson + RK2 | 34 | 33 | 30 | 27 | 34 | 158 | 4 | 1 |
| richardson + RK4 | 32 | 30 | 27 | 24 | 32 | 145 | 10 | 9 |
| OM + richardson | 32 | 33 | 34 | 23 | 34 | 156 | 8 | 0 |
| BE-chord | 19 | 18 | 17 | 14 | 21 | 89 | 24 | 47 |
| RK4 | 2 | 3 | 4 | 2 | 2 | 13 | 1 | 114 |
| scaled homotopy | 28 | 27 | 24 | 18 | 32 | 129 | 6 | 1 |
| newton homotopy | 15 | 15 | 15 | 15 | 15 | 75 | 10 | 56 |
