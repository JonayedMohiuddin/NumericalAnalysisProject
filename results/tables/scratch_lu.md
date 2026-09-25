**BE path + NR / FDXB with SciPy's sparse LU and with our dense LU**

| case | refiner | iterations (SuperLU) | iterations (our LU) | max \|x difference\| | time SuperLU (s) | time our LU (s) |
|---|---|---|---|---|---|---|
| tutorial (Sec. 3.3) | NR | 9 | 9 | 0.0e+00 | 0.014 | 0.007 |
| case69limit | NR | 1 | 1 | 3.6e-13 | 0.033 | 0.103 |
| case69limit | FDXB | 8 | 8 | 1.6e-13 | 0.029 | 0.095 |
