**BE path + NR / FDXB with SciPy's sparse LU and with our dense LU**

| case | refiner | iterations (SuperLU) | iterations (our LU) | max \|x difference\| | time SuperLU (s) | time our LU (s) |
|---|---|---|---|---|---|---|
| tutorial (Sec. 3.3) | NR | 9 | 9 | 0.0e+00 | 0.002 | 0.002 |
| case69limit | NR | 1 | 1 | 3.6e-13 | 0.004 | 0.011 |
| case69limit | FDXB | 8 | 8 | 1.6e-13 | 0.004 | 0.010 |
| case141limit | NR | 1 | 1 | 9.7e-14 | 0.005 | 0.039 |
| case141limit | FDXB | 8 | 8 | 4.7e-12 | 0.004 | 0.032 |
| case_ACTIVSg500limit | NR | 2 | 2 | 9.8e-15 | 0.010 | 3.454 |
| case_ACTIVSg500limit | FDXB | 13 | 13 | 1.2e-14 | 0.008 | 1.710 |
| case_ACTIVSg2000limit | NR | 3 | 3 | 4.3e-14 | 0.040 | 424.123 |
| case_ACTIVSg2000limit | FDXB | 22 | 22 | 6.2e-14 | 0.029 | 100.787 |
