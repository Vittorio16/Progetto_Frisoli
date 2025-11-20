import numpy as np
import sympy as sp

NUMERO_PENDOLI = 2

MASSE = np.array([1,1])
LUNGHEZZE = np.array([1,1])

state_vector = np.array([5,10,2,7])


# Describe all the symbols for simpy and lambdify
t = sp.symbols('t')

theta = [sp.Function(f"theta_{i+1}")(t) for i in range(NUMERO_PENDOLI)]
theta_dot = [sp.diff(theta[i],t) for i in range(NUMERO_PENDOLI)]
theta_dot_dot = [sp.diff(theta[i],(t,2)) for i in range(NUMERO_PENDOLI)]

thetaddot_sym = sp.symbols(f'thetaddot0:{NUMERO_PENDOLI}')


l = sp.symbols(f"l1:l{NUMERO_PENDOLI + 1}")
m = sp.symbols(f"m1:m{NUMERO_PENDOLI + 1}")

g,A,w = sp.symbols("g A w")


# Calculate the positions of each joint (x[-1] gets the last element)
# Use Sympy to simbolically write equations and solve them
x0 = A * sp.sin(w * t)

x = [x0 + l[0] * sp.sin(theta[0])]
y = [-l[0] * sp.cos(theta[0])]

for i in range(1,NUMERO_PENDOLI):
    x.append(x[-1] + l[i] * sp.sin(theta[i]))
    y.append(y[-1] - l[i] * sp.cos(theta[i]))
    
x_dot = [sp.diff(x[i], t) for i in range(NUMERO_PENDOLI)]
y_dot = [sp.diff(y[i], t) for i in range(NUMERO_PENDOLI)]


# Calculate energies and Lagrangian
T = sum(sp.Rational(1,2) * m[i] * (x_dot[i] ** 2 + y_dot[i] ** 2) for i in range(NUMERO_PENDOLI))
U = sum(m[i] * g * y[i] for i in range(NUMERO_PENDOLI))
L = T - U


eq = []
# Calulate Euler-Lagrange equations
for i in range(NUMERO_PENDOLI):
    dL_dtheta = sp.diff(L, theta[i])
    dL_dtheta_dot = sp.diff(L, theta_dot[i])
    dt_dL_dtheta_dot = sp.diff(dL_dtheta_dot, t)

    eq.append(sp.simplify(dt_dL_dtheta_dot - dL_dtheta))

eq_subs = [eq_i.subs({theta_dot_dot[j]: thetaddot_sym[j] for j in range(NUMERO_PENDOLI)}) for eq_i in eq]

# Find M and f, and the dd_theta as M^(-1)*f by solving the Euler-Lagrange equations
M = sp.zeros(NUMERO_PENDOLI)
f = sp.zeros(NUMERO_PENDOLI,1)

for i in range(NUMERO_PENDOLI):
    for j in range(NUMERO_PENDOLI):
        M[i,j] = eq[i].coeff(theta_dot_dot[j])
    f[i] = sp.simplify(eq[i] - sum(M[i,j] * theta_dot_dot[j] for j in range(NUMERO_PENDOLI)))
    
M_inv = sp.simplify(M.inv())
dd_theta = sp.simplify(M_inv * (-f))

ddtheta_fun = sp.lambdify(theta + theta_dot + [*l, *m, g], dd_theta, 'numpy')

# Suppose n = 3 pendulum segments
# Example numerical values:

import numpy as np

# Current angles (radians)
theta_vals = [0.1, 0.2]

# Current angular velocities
thetadot_vals = [0.0, 0.1]

# Lengths of each pendulum segment
l_vals = [1.0, 1.0]

# Masses of each pendulum segment
m_vals = [1.0, 1.0]

# Gravity
g_val = 9.81

# Combine all inputs in the order:
# theta_0, theta_1, theta_2, thetadot_0, thetadot_1, thetadot_2, l0, l1, l2, m0, m1, m2, g
args = theta_vals + thetadot_vals + l_vals + m_vals + [g_val]

# Call the function:
ddtheta_vals = ddtheta_fun(*args)

print(ddtheta_vals)
# for i in state_vector:
#     i = math.radians(i)
    
    
# Mtot = sum(MASSE)

# M11 = Mtot * (LUNGHEZZE[0] ** 2)
# M12 = MASSE[1] * LUNGHEZZE[0] * LUNGHEZZE[1] * math.cos(state_vector[0] - state_vector[1])
# M21 = MASSE[1] * LUNGHEZZE[0] * LUNGHEZZE[1] * math.cos(state_vector[0] - state_vector[1])
# M22 = MASSE[1] * (LUNGHEZZE[1] ** 2)

# mass_matrix = np.array([[M11,M12],[M21,M22]])