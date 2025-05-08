from sympy import *

xi, yi, ym, L, th_f, th_r, R_b, v_f, v_r = symbols('xi yi ym L th_f th_r R_b v_f v_r')
pi = 3.14159


# Define intermediate expressions
print(sqrt(1/3))
R_b = (L/sin(th_r+th_f))*sqrt(cos(th_r)*cos(th_f)*cos(th_f+th_r)+(1/4)*sin(th_r+th_f)**2)
psidot = ((v_f/cos(th_r))+(v_r/cos(th_f)))*sin(th_r+th_f)/(2*L)
v_b = R_b*psidot
pprint(simplify(v_b))
