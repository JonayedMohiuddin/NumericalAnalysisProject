**Setting S1: LU factorizations per case, or the outcome if not solved**

| case | paper | adaptive | PC | PC + adaptive | OM | OM + adaptive | OM + PC | OM + PC + adaptive | richardson | richardson + FE | richardson + RK2 | richardson + RK4 | OM + richardson | BE-chord | RK4 | scaled homotopy | newton homotopy |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| case18482 | 6 | 6 | 6 | 6 | 5 | 5 | 6 | 6 | 7 | 43 | 46 | 88 | 7 | other root | failed | 6 | failed |
| case27318 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 7 | 29 | 46 | 88 | 7 | other root | failed | 6 | 8 |
| case36964 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 13 | 46 | 46 | 88 | 13 | failed | failed | 6 | failed |
| case54636 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 7 | 29 | 46 | 88 | 7 | other root | failed | 6 | 8 |
| case109272 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 7 | 29 | 46 | 88 | 7 | other root | failed | 6 | other root |
| case69limit | 3 | 3 | 4 | 4 | 3 | 3 | 4 | 4 | 5 | 8 | 9 | 15 | 5 | 3 | 7 | 3 | 3 |
| case141limit | 3 | 3 | 4 | 4 | 3 | 3 | 4 | 4 | 5 | 8 | 9 | 15 | 5 | 3 | 7 | 3 | 3 |
| case_ACTIVSg500limit | 4 | 4 | 5 | 5 | 4 | 4 | 5 | 5 | 6 | 28 | 34 | 63 | 6 | 5 | failed | 4 | 5 |
| case_ACTIVSg2000limit | 5 | 14 | 6 | 6 | 5 | 14 | 6 | 6 | 6 | 34 | 51 | 87 | 6 | 6 | failed | 5 | 6 |
| case6024 | failed | failed | failed | failed | failed | failed | failed | failed | 73 | 106 | 112 | failed | 73 | failed | failed | failed | failed |
| case6243 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 7 | 43 | 46 | 88 | 7 | failed | failed | 6 | failed |
| case6748 | failed | 19 | failed | failed | failed | 18 | failed | failed | 49 | 52 | 52 | 100 | 49 | failed | failed | failed | failed |
| case7092 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 7 | 40 | 46 | 76 | 7 | 7 | failed | 6 | failed |
| case9961 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 7 | 40 | 46 | 88 | 7 | 7 | failed | 6 | failed |
| case10595 | 6 | 6 | 6 | 6 | 5 | 5 | 6 | 6 | 7 | 46 | 46 | 88 | 7 | other root | failed | 6 | failed |
| case12110 | 5 | 5 | 6 | 6 | 5 | 5 | 6 | 6 | 7 | 46 | 52 | 88 | 7 | 6 | failed | 6 | failed |
| case14limit | failed | failed | 13 | 13 | 10 | 10 | 10 | 10 | 14 | 34 | 57 | 129 | 12 | failed | failed | failed | failed |
| case_ieee30limit | 10 | 10 | 11 | 11 | 9 | 9 | 9 | 9 | 11 | 32 | 61 | 102 | 10 | 11 | failed | 10 | 11 |
| case57limit | 11 | 11 | 11 | 11 | 9 | 9 | 9 | 9 | 12 | 31 | 53 | 139 | 10 | 11 | failed | 10 | 11 |
| case89pegaselimit | 11 | 11 | 11 | 11 | 9 | 9 | 9 | 9 | 12 | 32 | 54 | 141 | 10 | 12 | failed | 11 | 12 |
| case118limit | 10 | 10 | 11 | 11 | 9 | 9 | 9 | 9 | 11 | 43 | 54 | other root | 9 | 11 | failed | 10 | 11 |
| case300limit | 10 | 10 | 11 | 11 | 8 | 8 | 9 | 9 | 11 | 37 | 49 | 78 | 9 | 11 | failed | 10 | 11 |
| case1354pegaselimit | 12 | 12 | 12 | 12 | 9 | 9 | 10 | 10 | 13 | 44 | 55 | 104 | 10 | 12 | failed | 12 | 12 |
| case2383wplimit | 10 | 10 | 11 | 19 | 9 | 9 | 10 | 18 | 11 | 40 | 54 | 78 | 10 | 11 | failed | 10 | 11 |
| case2736splimit | failed | failed | 13 | 13 | 11 | 11 | 9 | 9 | 14 | 45 | 51 | 105 | 12 | failed | failed | failed | failed |
| case2737soplimit | 12 | 12 | 12 | 12 | 10 | 10 | 10 | 10 | 13 | 47 | 62 | 92 | 11 | failed | failed | 12 | failed |
| case2746woplimit | 12 | 12 | 12 | 12 | 10 | 10 | 10 | 10 | 13 | 47 | 50 | 92 | 11 | 12 | failed | 12 | failed |
| case2746wplimit | 12 | 12 | 12 | 12 | 10 | 10 | 10 | 10 | 13 | 44 | 56 | 116 | 11 | failed | failed | 12 | failed |
| case2869pegaselimit | 10 | 10 | 11 | 11 | 9 | 9 | 10 | 10 | 11 | 42 | 54 | 90 | 10 | 11 | failed | 10 | 11 |
| case3012wplimit | other root | 19 | 11 | 11 | other root | 18 | 10 | 10 | other root | 53 | 56 | 116 | other root | failed | failed | other root | failed |
| case3120splimit | 12 | 12 | 12 | 12 | 10 | 10 | 10 | 10 | 13 | 47 | 56 | 92 | 10 | failed | failed | 12 | failed |
| case3375wplimit | other root | other root | 11 | 11 | other root | other root | 9 | 9 | other root | 53 | 50 | 91 | other root | 12 | failed | other root | failed |
| case9241pegaselimit | 11 | 11 | 11 | 11 | 9 | 9 | 10 | 10 | 46 | 52 | 55 | 91 | 45 | 11 | failed | 11 | 12 |
| case13659pegaselimit | 8 | 8 | 9 | 9 | 8 | 8 | 8 | 8 | 9 | 33 | 48 | 90 | 9 | 10 | failed | 8 | other root |
