import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint
import sympy as sp

def spring(start, end, nodes, width):
    """!
    Return a list of points corresponding to a spring.
    @param r1 (array-like) The (x, y) coordinates of the first endpoint.
    @param r2 (array-like) The (x, y) coordinates of the second endpoint.
    @param nodes (int) The number of spring "nodes" or coils.
    @param width (int or float) The diameter of the spring.
    @return An array of x coordinates and an array of y coordinates.
    """

    # Check that nodes is at least 1.
    nodes = max(int(nodes), 1)

    # Convert to numpy array to account for inputs of different types/shapes.
    start, end = np.array(start).reshape((2,)), np.array(end).reshape((2,))

    # If both points are coincident, return the x and y coords of one of them.
    if (start == end).all():
        return start[0], start[1]

    # Calculate length of spring (distance between endpoints).
    length = np.linalg.norm(np.subtract(end, start))

    # Calculate unit vectors tangent (u_t) and normal (u_t) to spring.
    u_t = np.subtract(end, start) / length
    u_n = np.array([[0, -1], [1, 0]]).dot(u_t)

    # Initialize array of x (row 0) and y (row 1) coords of the nodes+2 points.
    spring_coords = np.zeros((2, nodes + 2))
    spring_coords[:,0], spring_coords[:,-1] = start, end

    # Check that length is not greater than the total length the spring
    # can extend (otherwise, math domain error will result), and compute the
    # normal distance from the centerline of the spring.
    normal_dist = math.sqrt(max(0, width**2 - (length**2 / nodes**2))) / 2

    # Compute the coordinates of each point (each node).
    for i in range(1, nodes + 1):
        spring_coords[:,i] = (
            start
            + ((length * (2 * i - 1) * u_t) / (2 * nodes))
            + (normal_dist * (-1)**i * u_n))

    return spring_coords[0,:], spring_coords[1,:]

def odesys(y, t, m1, m, c, g, r, c1):
    dy = np.zeros(4)
    dy[0] = y[2]
    dy[1] = y[3]

    a11 = (2 * m1 + m) * r ** 2 + m * y[1] ** 2 + 2 * m * r * y[1] * sp.sin(y[0])
    a12 = -m * r * sp.cos(y[0])
    a21 = -r * sp.cos(y[0])
    a22 = 1

    b1 = -c1 * r ** 2 * y[0] - m * g * y[1] * sp.cos(y[0]) - m * r * sp.cos(y[0]) * y[2] ** 2 * y[1] - 2 * m * (
            y[1] + r * sp.sin(y[0])) * y[3] * y[2]
    b2 = -2 * (c / m) * y[1] - g * sp.sin(y[0]) + y[1] * y[2] ** 2

    dy[2] = (b1 * a22 - b2 * a12) / (a11 * a22 - a12 * a21)
    dy[3] = (b2 * a11 - b1 * a21) / (a11 * a22 - a12 * a21)

    return dy


# data of task
m1 = 1
m = 0.5
r = 1
c = 5
c1 = 5
t0 = 0
g = 9.81

t_fin = 20

T = np.linspace(0, t_fin, 1001)

phi0 = math.pi / 2
s0 = 0
dphi0 = 1
ds0 = 0

y0 = [phi0, s0, dphi0, ds0]

Y = odeint(odesys, y0, T, (m1, m, c, g, r, c1))

print(Y.shape)

phi = Y[:, 0]
s = Y[:, 1]
# Ввод переменной t и радиусов необходимых окружностей + ввод угла поворота шариков
t = sp.Symbol('t')

# Построение графика и подграфика с выравниванием осей
fig = plt.figure(figsize=(17, 4))
ax1 = fig.add_subplot(1, 2, 1)
ax1.set(xlim=[-15, 15], ylim=[-5, 5])
ax1.axis('equal')

R = r
L = 0.5

# points
X_O = L + R + R * np.cos(phi)
Y_O = R

line1 = ax1.plot([0, 0], [0, 2 * R], color='black')[0]
line2 = ax1.plot([0, 2 * R + L], [0, 0], color='black')[0]
circle1 = plt.Circle((X_O[0], Y_O), R, color='black', fill=False)
cargo = plt.Rectangle((X_O[0] - s[0] * (math.sqrt(2) / 2), Y_O + s[0] * (math.sqrt(2) / 2)), width=0.2, height=0.2, angle=45)

line3 = ax1.plot([X_O[0] - R * (math.sqrt(2) / 2), X_O[0] + R * (math.sqrt(2) / 2)],
                 [Y_O + R * (math.sqrt(2) / 2), Y_O - R * (math.sqrt(2) / 2)], color='black', linewidth=13, alpha=0.1)[
    0]

# spring1
Sh = 0.5
b = 1 / (c1 - 2)
X_Spr = np.zeros(c1)
Y_Spr = np.zeros(c1)
X_Spr[0] = 0
Y_Spr[0] = 0
X_Spr[c1 - 1] = 1
Y_Spr[c1 - 1] = 0
for i in range(c1 - 2):
    X_Spr[i + 1] = b * (i + 1) - b / 2
    Y_Spr[i + 1] = Sh * (-1) ** i

Y_G = R
L_Spr = X_O
spring1 = ax1.plot(X_Spr * L_Spr[0], Y_Spr + Y_G)[0]

# spring2
spring2 = ax1.plot(*spring((X_O[0] - R * (math.sqrt(2) / 2), Y_O + R * (math.sqrt(2) / 2)), ((X_O[0] - s[0] * (math.sqrt(2) / 2), Y_O + s[0] * (math.sqrt(2) / 2))), 12, 0.2))[0]


#spring3
spring3 = ax1.plot(*spring((X_O[0] + R * (math.sqrt(2) / 2), Y_O - R * (math.sqrt(2) / 2)), ((X_O[0] - s[0] * (math.sqrt(2) / 2), Y_O + s[0] * (math.sqrt(2) / 2))), 12, 0.2))[0]

ax1.add_patch(circle1)
ax1.add_patch(cargo)


def anime(i):
    circle1.center = X_O[i], Y_O
    cargo.set_xy((X_O[i] - s[i] * (math.sqrt(2) / 2), Y_O + s[i] * (math.sqrt(2) / 2)))
    spring1.set_data((X_Spr * L_Spr[i], Y_Spr + Y_G))
    spring2.set_data(*spring((X_O[i] - R * (math.sqrt(2) / 2), Y_O + R * (math.sqrt(2) / 2)), (
    X_O[i] - s[i] * (math.sqrt(2) / 2), Y_O + s[i] * (math.sqrt(2) / 2)), 12, 0.2))
    spring3.set_data(*spring((X_O[i] + R * (math.sqrt(2) / 2), Y_O - R * (math.sqrt(2) / 2)),
                             (X_O[i] - s[i] * (math.sqrt(2) / 2), Y_O + s[i] * (math.sqrt(2) / 2)), 12, 0.2))
    line3.set_data([X_O[i] - R * (math.sqrt(2) / 2), X_O[i] + R * (math.sqrt(2) / 2)],
                   [Y_O + R * (math.sqrt(2) / 2),
                    Y_O - R * (math.sqrt(2) / 2)])
    return circle1, spring1, line3, spring2, spring3


anim = FuncAnimation(fig, anime, frames=500, interval=10)
plt.show()
