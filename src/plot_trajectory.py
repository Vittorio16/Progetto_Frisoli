import numpy as np
from calculate_trajectory import calc_trajectory, A_val, w_val, NUMERO_PENDOLI, LUNGHEZZE
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Computes position with calc_trajectory, then transforms it in arrays x = [x_i,t_i], y = [y_i, t_i]
def compute_positions():
    theta_array = solution.y[:NUMERO_PENDOLI]
    
    x = np.zeros((NUMERO_PENDOLI + 1, theta_array.shape[1]))
    y = np.zeros((NUMERO_PENDOLI + 1, theta_array.shape[1]))
    
    # Position of the first joint - solution only returns angles, so need to recalculate this
    x[0] = A_val * np.sin(w_val * solution.t)
    y[0] = 0
    
    # First pendulum
    x[1] = x[0] + LUNGHEZZE[0] * np.sin(theta_array[0])
    y[1] = y[0] - LUNGHEZZE[0] * np.cos(theta_array[0])

    # Chain the links
    for i in range(2, NUMERO_PENDOLI + 1):
        x[i] = x[i-1] + LUNGHEZZE[i-1] * np.sin(theta_array[i-1])
        y[i] = y[i-1] - LUNGHEZZE[i-1] * np.cos(theta_array[i-1])
    
    return x, y
    return x, y
    
    
# Creates the actual animation using matplotlib
def update(frame):
    xs = list(x[:, frame])
    ys = list(y[:, frame])
    line.set_data(xs, ys)
    return line,



t_span = (0, 20)
t_eval = np.linspace(t_span[0], t_span[1], 5000)

solution = calc_trajectory(t_span, t_eval)

x, y = compute_positions()
print(y[0,:])

fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlim(-sum(LUNGHEZZE)-A_val-0.5, sum(LUNGHEZZE)+A_val+0.5)
ax.set_ylim(-sum(LUNGHEZZE)-0.5, 1)

line, = ax.plot([], [], 'o-', lw=2)

ani = FuncAnimation(fig, update, frames=len(solution.t), interval=10, blit=True)
plt.show()