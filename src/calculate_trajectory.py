import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

# Check units
g_val = 9.81
A_val = 2
w_val = 2 * np.pi * 0.4

initial_time = 0
initial_state_vector = np.radians(np.array([90,0,0,0]))

NUMERO_PENDOLI = 2
MASSE = np.array([1000,1])
LUNGHEZZE = np.array([1,0.2])

# Describe all the symbols for simpy and lambdify
t = sp.symbols('t')

theta = [sp.Function(f"theta_{i+1}")(t) for i in range(NUMERO_PENDOLI)]
theta_dot = [sp.diff(theta[i],t) for i in range(NUMERO_PENDOLI)]

# theta_dot_dot = [sp.diff(theta[i],(t,2)) for i in range(NUMERO_PENDOLI)]

sym_theta_dot_dot = sp.symbols(f"theta_dot_dot_1:{NUMERO_PENDOLI + 1}")
l = sp.symbols(f"l1:{NUMERO_PENDOLI + 1}")
m = sp.symbols(f"m1:{NUMERO_PENDOLI + 1}")
g,A,w = sp.symbols("g A w")


# Calculate the positions of each joint (x[-1] gets the last element) - each saved in the array except for the first one
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

M = sp.simplify(M)
f = sp.simplify(f)

# M_inv = sp.simplify(M.inv())
# dd_theta = sp.simplify(M_inv * (f))

# Translates from symbolic to numerical, passing t for the forced joint position, and all the other necessary elements
# It calculates dd_theta at time t

# ddtheta_fun = sp.lambdify([t] + theta + theta_dot + [*l, *m, g, A, w], dd_theta, 'numpy')
M_fun = sp.lambdify(theta + [*l, *m], M, "numpy")
f_fun = sp.lambdify([t] + theta + theta_dot + [*l, *m, g, A, w], f, "numpy")


# Given time and state, returns the theta_dot and theta_dot_dot array
def pendulum_ode(t, state):
    theta_vals = state[:NUMERO_PENDOLI]
    theta_dot_vals = state[NUMERO_PENDOLI:]

    # args = [t] + list(theta_vals) + list(theta_dot_vals) \
    #            + list(LUNGHEZZE) + list(MASSE) + [g_val, A_val, w_val]

    # dd_theta_vals = np.array(ddtheta_fun(*args)).flatten()
    M_num = M_fun(*theta_vals, *LUNGHEZZE, *MASSE)
    f_num = f_fun(t, *theta_vals, *theta_dot_vals, *LUNGHEZZE, *MASSE, g_val, A_val, w_val)

    dd_theta_vals = np.linalg.solve(M_num, f_num).flatten()

    # return state derivative
    return np.concatenate([theta_dot_vals, dd_theta_vals])


# It integrates over the time span, returning samples according to t_eval
def calc_trajectory(t_span, t_eval):
    
    solution = solve_ivp(
        pendulum_ode,
        t_span,
        initial_state_vector,
        t_eval=t_eval,
        method='RK45',
        rtol=1e-9,
        atol=1e-9
    )
    
    return solution
