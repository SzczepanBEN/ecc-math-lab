# Experiment A summary (auto-generated)

Permutations per null: 200; seed 1.

## Ordering statistics (half period, Monte-Carlo p-values vs random permutation)

| curve | acf_max_abs | ljung_box_Q | diff_chi2 | serial2d_chi2 | rank_corr_k_abs | ks_first_1000 | periodogram_peak_ratio |
|---|---|---|---|---|---|---|---|
| GEN-A | 0.776 | 0.622 | 0.846 | 0.537 | 0.910 | 0.687 | 0.100 |
| GEN-B | 0.398 | 0.453 | 0.985 | 0.692 | 0.080 | 0.756 | 0.483 |
| GEN-C6 | 0.960 | 0.965 | 0.930 | 0.189 | 0.572 | 0.254 | 0.144 |
| J0-A | 0.657 | 0.547 | 0.761 | 0.453 | 0.413 | 0.388 | 0.408 |
| J0-B | 0.821 | 0.468 | 0.015 | 0.677 | 0.547 | 0.985 | 0.303 |
| J0-SS | 0.701 | 0.706 | 0.781 | 0.104 | 0.642 | 0.239 | 0.299 |
| J1728 | 0.224 | 0.348 | 0.219 | 0.468 | 0.975 | 0.861 | 0.781 |

1 of 49 p-values < 0.05 (expected under the null: 2.5); min p = 0.0149; Bonferroni threshold 0.00102.

## Value distribution and exact symmetries

Spread p-values are (upper/lower) tails; 'naive' = arbitrary relabelling of F_p, 'sym' = relabelling commuting with v -> -v (correct null when the occupancy is negation-invariant).

| curve | x occurrence profile | y occurrence profile | spread-x naive | spread-x sym | spread-y naive | spread-y sym | x_{n-k}=x_k | y_{n-k}=-y_k |
|---|---|---|---|---|---|---|---|---|
| GEN-A | {0: 64110, 2: 64373} | {0: 42601, 1: 64450, 3: 21432} | 0.612/0.393 | n/a | 0.025/0.980 | 0.065/0.940 | 1.000 | 1.000 |
| GEN-B | {0: 49613, 2: 49764} | {0: 33059, 1: 49712, 2: 2, 3: 16604} | 0.010/0.995 | n/a | 0.612/0.393 | 0.567/0.438 | 1.000 | 1.000 |
| GEN-C6 | {0: 61964, 1: 1, 2: 61922} | {0: 41342, 1: 61895, 3: 20650} | 0.537/0.468 | n/a | 0.995/0.010 | 0.945/0.060 | 1.000 | 1.000 |
| J0-A | {0: 58870, 2: 59553} | {0: 78721, 3: 39702} | 0.308/0.697 | n/a | 0.841/0.164 | 0.632/0.373 | 1.000 | 1.000 |
| J0-B | {0: 56998, 2: 57549} | {0: 76181, 3: 38366} | 0.453/0.552 | n/a | 0.055/0.950 | 0.080/0.925 | 1.000 | 1.000 |
| J0-SS | {0: 36134, 1: 1, 2: 36134} | {1: 72269} | 0.428/0.577 | n/a | 1.000/1.000 | 1.000/1.000 | 1.000 | 1.000 |
| J1728 | {0: 35896, 1: 1, 2: 36372} | {0: 23848, 1: 36257, 2: 4, 3: 12160} | 0.249/0.756 | 0.353/0.652 | 0.851/0.154 | 0.761/0.244 | 1.000 | 1.000 |

## Group-structure and automorphism checks

- **GEN-A** torsion: {}; automorphisms: None
  - y-profile vs Chebotarev: {'observed_0_1_3': [42601, 64450, 21432], 'expected': [42827.7, 64241.5, 21413.8], 'chi2_p_value': 0.38833919469730305, 'count_mult_2': 0}
- **GEN-B** torsion: {}; automorphisms: None
  - y-profile vs Chebotarev: {'observed_0_1_3': [33059, 49712, 16604], 'expected': [33125.0, 49687.5, 16562.5], 'chi2_p_value': 0.8835709827321102, 'count_mult_2': 2}
- **GEN-C6** torsion: {'k=n/2: y==0': True, '#k with y_k=0': 1, 'k=n/3,2n/3: psi_3(x)==0': True, '#k with psi_3(x_k)=0': 1}; automorphisms: None
  - y-profile vs Chebotarev: {'observed_0_1_3': [41342, 61895, 20650], 'expected': [41295.7, 61943.5, 20647.8], 'chi2_p_value': 0.9559080881678499, 'count_mult_2': 0}
- **J0-A** torsion: {}; automorphisms: [{'map': 'x->beta x, y->y', 'lambda': 4868, 'lambda_minimal_poly_check': {'l^2+l+1 mod n': 0, 'l^2+1 mod n': 114239}, 'fraction_of_k_satisfied': 1.0}, {'map': 'x->beta x, y->y', 'lambda': 114238, 'lambda_minimal_poly_check': {'l^2+l+1 mod n': 0, 'l^2+1 mod n': 4869}, 'fraction_of_k_satisfied': 1.0}]
  - x^3 occurrence profile: {0: 98572, 6: 19851}
- **J0-B** torsion: {}; automorphisms: [{'map': 'x->beta x, y->y', 'lambda': 95325, 'lambda_minimal_poly_check': {'l^2+l+1 mod n': 0, 'l^2+1 mod n': 19774}, 'fraction_of_k_satisfied': 1.0}, {'map': 'x->beta x, y->y', 'lambda': 19773, 'lambda_minimal_poly_check': {'l^2+l+1 mod n': 0, 'l^2+1 mod n': 95326}, 'fraction_of_k_satisfied': 1.0}]
  - x^3 occurrence profile: {0: 95364, 6: 19183}
- **J0-SS** torsion: {'k=n/2: y==0': True, '#k with y_k=0': 1, 'k=n/3,2n/3: psi_3(x)==0': True, '#k with psi_3(x_k)=0': 1}; automorphisms: None
  - x^3 occurrence profile: {0: 36134, 1: 1, 2: 36134}
- **J1728** torsion: {'k=n/2: y==0': True, '#k with y_k=0': 1}; automorphisms: [{'map': 'x->-x, y->i y', 'lambda': 19785, 'lambda_minimal_poly_check': {'l^2+l+1 mod n': 19785, 'l^2+1 mod n': 0}, 'fraction_of_k_satisfied': 1.0}]
