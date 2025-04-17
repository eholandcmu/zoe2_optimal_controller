from sympy import *

xi, yi, ym, L, th_f, th_r, R_b = symbols('xi yi ym L th_f th_r R_b')
pi = 3.14159


# Define intermediate expressions
yi_expr = (L * tan(th_r)) / (tan(th_r) + tan(th_f))
xi_expr = (-L) / (tan(th_r) + tan(th_f))
ym_expr = L / 2

# Define the final equation using those expressions
R_b_expr = sqrt((yi_expr - ym_expr)**2 + xi_expr**2)

# Optional: simplify
R_b_simplified = simplify(R_b_expr)

pprint(R_b_simplified)
