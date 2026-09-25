**Choosing delta: rules fitted on the paper's 9 cases, tested on the other 25. The golden search row excludes the cost of the search itself**

| how delta is chosen | solved, training | mean LUs, training | solved, test | mean LUs, test |
|---|---|---|---|---|
| paper, delta = 0.02 | 9 of 9 | 5.00 | 19 of 25 | 9.47 |
| fixed rule, delta = 0.0111 | 9 of 9 | 5.11 | 19 of 25 | 9.53 |
| spectral rule, delta = 1.86 lambda_min | 5 of 9 | 5.20 | 20 of 25 | 11.00 |
| golden search per case | 9 of 9 | 4.89 | 23 of 25 | 9.13 |
