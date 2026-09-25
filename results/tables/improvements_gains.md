**Cases the paper's method does not solve that another configuration solves**

| setting | case | paper's method with 30 NR iterations | solved within 10 by |
|---|---|---|---|
| S1 paper | case6024 | diverges | richardson; richardson + FE; richardson + RK2; OM + richardson |
| S1 paper | case6748 | diverges | adaptive; OM + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S1 paper | case14limit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S1 paper | case2736splimit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S1 paper | case3012wplimit | reaches another root | adaptive; PC; PC + adaptive; OM + adaptive; OM + PC; OM + PC + adaptive; richardson + FE; richardson + RK2; richardson + RK4 |
| S1 paper | case3375wplimit | reaches another root | PC; PC + adaptive; OM + PC; OM + PC + adaptive; richardson + FE; richardson + RK2; richardson + RK4; BE-chord |
| S2 Sec 4.4 | case36964 | diverges | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S2 Sec 4.4 | case6024 | diverges | adaptive; OM + adaptive; richardson; richardson + FE; OM + richardson |
| S2 Sec 4.4 | case6748 | diverges | adaptive; OM + adaptive; richardson; richardson + RK2 |
| S2 Sec 4.4 | case14limit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S2 Sec 4.4 | case2736splimit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S2 Sec 4.4 | case3012wplimit | reaches another root | adaptive; PC; PC + adaptive; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson; BE-chord |
| S2 Sec 4.4 | case3375wplimit | reaches another root | PC; PC + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson; BE-chord |
| S3 Sec 4.4 | case36964 | diverges | PC; PC + adaptive; OM + PC; OM + PC + adaptive; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S3 Sec 4.4 | case54636 | reaches another root | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; OM + richardson |
| S3 Sec 4.4 | case109272 | reaches another root | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; OM + richardson |
| S3 Sec 4.4 | case6024 | diverges | richardson; richardson + FE; OM + richardson |
| S3 Sec 4.4 | case6748 | diverges | richardson + RK2; OM + richardson |
| S3 Sec 4.4 | case14limit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + RK2; richardson + RK4; OM + richardson |
| S3 Sec 4.4 | case2736splimit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S3 Sec 4.4 | case3012wplimit | reaches another root | adaptive; PC; PC + adaptive; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson; BE-chord |
| S3 Sec 4.4 | case3375wplimit | reaches another root | PC; PC + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson; BE-chord |
| S4 weak K | case27318 | reaches another root | PC; PC + adaptive; OM + PC; OM + PC + adaptive; richardson + FE; richardson + RK2; richardson + RK4; BE-chord; newton homotopy |
| S4 weak K | case54636 | reaches another root | PC; PC + adaptive; OM + PC; OM + PC + adaptive; richardson + FE; richardson + RK2; richardson + RK4; newton homotopy |
| S4 weak K | case109272 | reaches another root | richardson + FE; richardson + RK2; richardson + RK4 |
| S4 weak K | case7092 | diverges | richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S4 weak K | case10595 | diverges | richardson; richardson + FE; richardson + RK2; OM + richardson |
| S4 weak K | case12110 | diverges | richardson; richardson + FE; richardson + RK2; OM + richardson |
| S4 weak K | case14limit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S4 weak K | case2736splimit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S4 weak K | case13659pegaselimit | reaches another root | richardson + FE; richardson + RK2; richardson + RK4; BE-chord |
| S5 strong K | case27318 | reaches another root | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson; BE-chord; scaled homotopy; newton homotopy |
| S5 strong K | case54636 | reaches another root | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson; BE-chord; scaled homotopy; newton homotopy |
| S5 strong K | case109272 | reaches another root | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson; scaled homotopy |
| S5 strong K | case14limit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
| S5 strong K | case2736splimit | slow (converges in 11) | PC; PC + adaptive; OM; OM + adaptive; OM + PC; OM + PC + adaptive; richardson; richardson + FE; richardson + RK2; richardson + RK4; OM + richardson |
