# Numbers copied from the paper so the scripts can print them next to ours.
# None stands for "-" in the paper.

# Tables 2-4: ||g|| at t = 0, 0.005, 0.01, 0.51, 1.0, then after NR iterations 1, 2, 3
TABLE2_FE = {
    "case18482": [532, 65, 9.7, 2e3, 2e3, None, None, None],
    "case27318": [200, 25, 1.2, 27, 2e7, None, None, None],
    "case36964": [532, 65, 3.5, 2e3, 3e5, None, None, None],
    "case54636": [200, 25, 1.2, 27, 1e9, None, None, None],
    "case109272": [200, 25, 1.2, 27, 1e9, None, None, None],
    "case69limit": [0.012, 2e-4, 1e-5, 5e-4, 2e-5, 2e-8, 3e-12, None],
    "case141limit": [0.006, 2e-4, 9e-6, 5e-4, 2e-5, 4e-8, 1e-11, None],
    "case_ACTIVSg500limit": [13, 1.2, 0.02, 0.9, 1.0, 8e-3, 2e-6, 4e-13],
    "case_ACTIVSg2000limit": [54, 3.9, 0.8, 210, 7e3, None, None, None],
}
TABLE3_RK2 = {
    "case18482": [532, 65, 33, 3e7, 1e11, None, None, None],
    "case27318": [200, 25, 13, 6e5, 1e10, None, None, None],
    "case36964": [532, 65, 33, 3e9, 2e14, None, None, None],
    "case54636": [200, 25, 13, 7e5, 1e10, None, None, None],
    "case109272": [200, 25, 13, 1e7, 5e12, None, None, None],
    "case69limit": [0.012, 2e-4, 1e-4, 9e-5, 5e-5, 1e-9, None, None],
    "case141limit": [0.006, 2e-4, 8e-5, 4e-4, 2e-4, 7e-9, None, None],
    "case_ACTIVSg500limit": [13, 1.2, 0.6, 1e4, 5e3, None, None, None],
    "case_ACTIVSg2000limit": [54, 3.9, 1.9, 1e4, 4e4, None, None, None],
}
TABLE4_BE = {
    "case18482": [532, 65, 33, 7.2, 3.7, 0.04, 1.5e-5, 1e-11],
    "case27318": [200, 25, 13, 0.5, 0.25, 0.1, 2.5e-5, 7e-12],
    "case36964": [532, 65, 33, 6.9, 3.6, 0.27, 6e-4, 4e-9],
    "case54636": [200, 25, 13, 0.5, 0.25, 0.1, 3e-5, 7e-12],
    "case109272": [200, 25, 13, 0.5, 0.25, 0.1, 3e-5, 7e-12],
    "case69limit": [0.012, 2e-4, 1e-4, 2e-6, 5e-7, 5e-12, None, None],
    "case141limit": [0.006, 2e-4, 8e-5, 1e-6, 5e-7, 6e-11, None, None],
    "case_ACTIVSg500limit": [13, 1.2, 0.6, 0.02, 9e-3, 1e-6, 2e-13, None],
    "case_ACTIVSg2000limit": [54, 3.9, 2, 2.8, 1.5, 0.01, 8e-7, 9e-13],
}

# Table 5: ||g|| at t = 0, 0.005, 0.01, 0.02, 0.12, 0.62, 1.0, then NR iterations 1, 2, 3
TABLE5 = {
    ("case36964", "FE"): [532, 65, 3.5, 0.7, 5.7, 57, 114, 28, 5, 0.2],
    ("case36964", "RK2"): [532, 65, 33, 16, 9.6, 4.2, 2.6, 0.2, 3e-4, 1e-9],
    ("case36964", "BE"): [532, 65, 33, 16, 4.0, 1.4, 0.9, 0.09, 9e-5, 1e-10],
    ("case109272", "FE"): [200, 25, 1.2, 0.02, 0.07, 0.3, 0.17, 0.03, 6e-4, 2e-9],
    ("case109272", "RK2"): [200, 25, 13, 6.3, 1.2, 0.5, 0.3, 0.19, 2e-4, 2e-9],
    ("case109272", "BE"): [200, 25, 13, 6.3, 1.1, 0.2, 0.1, 0.02, 1e-6, 1e-11],
    ("case_ACTIVSg2000limit", "FE"): [55, 3.9, 0.8, 0.7, 3.6, 58, 91, 16, 1.0, 6e-3],
    ("case_ACTIVSg2000limit", "RK2"): [55, 3.9, 1.9, 1.0, 5.0, 5.5, 3.4, 0.2, 3e-4, 1e-9],
    ("case_ACTIVSg2000limit", "BE"): [55, 3.9, 1.9, 1.0, 1.0, 0.3, 0.2, 5e-4, 3e-9, None],
}

# Table 6: CPU time in seconds for case109272, case36964, case_ACTIVSg2000limit
TABLE6 = {
    ("NR-MAT", None): [3.8889, 1.8862, 0.3100],
    ("FE", 0.20): [5.6750, 3.3580, 0.2120],
    ("FE", 0.25): [5.6500, "fail", 0.2340],
    ("FE", 0.30): [6.6210, 3.3680, 0.2240],
    ("RK2", 0.20): [8.9580, 4.3800, 0.2550],
    ("RK2", 0.25): [8.7950, 3.5900, 0.2610],
    ("RK2", 0.30): [9.2150, 3.6770, 0.2270],
    ("BE", 0.20): [6.6550, 2.7160, 0.1890],
    ("BE", 0.25): [6.6540, 2.7110, 0.2500],
    ("BE", 0.30): [6.5690, 2.7860, 0.2900],
}

# Table 7: NR-MAT, NR-flat, dh1, delta, GSH iter, BE(NR), RK2(NR), BE(FDXB), RK2(FDXB)
TABLE7 = {
    "case18482": [6, "fail", 0.5, 1, 10, 4, 3, 33, 29],
    "case27318": [5, "fail", 0.25, 0.125, 25, 4, 3, 29, 25],
    "case36964": [6, "fail", 0.25, 1, 27, 4, 3, 33, 29],
    "case54636": [5, "fail", 0.01, 0.01, "fail", 4, 3, 29, 25],
    "case109272": [5, "fail", 0.01, 0.01, "fail", 4, 3, 19, 25],
    "case69limit": [4, 4, 1, 1, 4, 1, 1, 8, 8],
    "case141limit": [3, 3, 1, 1, 3, 1, 1, 8, 6],
    "case_ACTIVSg500limit": [3, 4, 1, 1, 4, 2, 2, 13, 11],
    "case_ACTIVSg2000limit": [3, 5, 1, 1, 5, 3, 3, 22, 20],
}

# Table 8: CPU time in seconds
TABLE8_CASES = ["case109272", "case54636", "case36964", "case27318", "case18482"]
TABLE8 = {
    "NR-MAT": [3.8889, 2.1236, 1.8862, 0.8912, 0.8241],
    "GSH-NR": ["fail", "fail", 40.6900, 21.5540, 9.4790],
    "BE(NR)": [4.2730, 2.1840, 1.7260, 1.1430, 0.9040],
    "RK2(NR)": [6.3920, 3.1070, 4.1100, 1.6240, 1.1630],
    "BE(FDXB)": [2.6890, 1.4070, 1.0840, 0.7560, 0.6250],
    "RK2(FDXB)": [4.1690, 2.0840, 1.5890, 1.0600, 0.9160],
}

# Figure 1a: NR iterations to reach 1e-5 from x0 = [1, 1 - eps]
FIG1A_ITERS = {0.005: 11, 0.01: 10, 0.05: 9}
# Figures 1(b)-(d) and 2, read off the plots: (dt, K) -> NR iterations after the homotopy
FIG1_NR_ITERS = {
    (0.5, 0.05): {"FE": 6, "BE": 3, "RK2": 6},
    (0.5, 0.01): {"FE": 4, "BE": 2, "RK2": 4},
    (0.5, 0.005): {"FE": 6, "BE": 3, "RK2": 7},
    (0.125, 0.05): {"FE": 3, "BE": 2, "RK2": 3},
}
