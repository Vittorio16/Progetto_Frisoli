import numpy as np
import sympy as sp

NUMERO_PENDOLI = 2
MASSE = np.array([1,1])
LUNGHEZZE = np.array([1,1])

# Describe all the symbols for simpy and lambdify
t = sp.symbols('t')

theta = [sp.Function(f"theta_{i+1}")(t) for i in range(NUMERO_PENDOLI)]
theta_dot = [sp.diff(theta[i],t) for i in range(NUMERO_PENDOLI)]

# theta_dot_dot = [sp.diff(theta[i],(t,2)) for i in range(NUMERO_PENDOLI)]

sym_theta_dot_dot = sp.symbols(f"theta_dot_dot_1:{NUMERO_PENDOLI + 1}")
l = sp.symbols(f"l1:{NUMERO_PENDOLI + 1}")
m = sp.symbols(f"m1:{NUMERO_PENDOLI + 1}")
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


# Write the equation symbolically
subs_dict = {sp.diff(theta[j], (t,2)): sym_theta_dot_dot[j] for j in range(NUMERO_PENDOLI)}
eq_subs = [eq_i.subs(subs_dict) for eq_i in eq]

# Find M and f, and the dd_theta as M^(-1)*f by solving the Euler-Lagrange equations    
M, f = sp.linear_eq_to_matrix(eq_subs, sym_theta_dot_dot)
f = -f  # move RHS to the form M * dd_theta = f

M = sp.simplify(M)
f = sp.simplify(f)

M_inv = sp.simplify(M.inv())
dd_theta = sp.simplify(M_inv * (-f))

# Translates from symbolic to numerical, passing t for the forced joint position, and all the other necessary elements
# It calculates dd_theta at time t
ddtheta_fun = sp.lambdify([t] + theta + theta_dot + [*l, *m, g, A, w], dd_theta, 'numpy')


def pendulum_ode(time, X, ddtheta_fun):
    
    g_val = 9.81
    A_val = 2
    w_val = np.radians(5)

    args = [time] + list(X) + list(LUNGHEZZE) + list(MASSE) + [g_val, A_val, w_val]
    
    dd_theta_vals = np.array(ddtheta_fun(*args)).flatten



initial_state_vector = np.radians(np.array([5,10,2,7]))
t_start = 0



