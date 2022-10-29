import math
import numpy as np
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def odesys(y, t, m1, m2, c, l, e, alpha, g):

    dy = np.zeros(4)
    dy[0] = y[2]
    dy[1] = y[3]

    a11 = ((5 / 6) * m1 + (4 / 3) * m2) * r * r
    a12 = 2 * m2 * l * e * math.sin(alpha)
    a21 = 2 * e * math.sin(alpha)
    a22 = l

    b1 = (m1 * np.sin(y[0]) + (2 * m2 * np.cos(y[0] - math.pi / 6))) * e * g - c * y[0] + 2 * m2 * l * e * y[3] ** 2 * math.cos(alpha)
    b2 = -g * np.sin(y[1]) - 2 * e * y[2] ** 2 * math.cos(alpha)

    dy[2] = (b1 * a22 - b2 * a12) / (a11 * a22 - a12 * a21)
    dy[3] = (b2 * a11 - b1 * a21) / (a11 * a22 - a12 * a21)

    return dy


# data of task
l = 5
alpha = 30
m1 = 2
m2 = 0
r = 5
c = 10
e = r / math.sqrt(3)
g = 9.81

t_fin = 20

T = np.linspace(0, t_fin, 1001)

phi0 = math.pi / 10
tau0 = 0.1
dphi0 = math.pi / 7
dtau0 = 0.3

y0 = [phi0, tau0, dphi0, dtau0]

Y = odeint(odesys, y0, T, (m1, m2, c, l, e, alpha, g))

print(Y.shape)

phi = Y[:, 0]
tau = Y[:, 1]

# CODE FROM OLD LAB

X_C = 0
Y_C = 0
X_O = math.sqrt(e / 2)
Y_O = - math.sqrt(e / 2)
X_A = np.sin(phi) * r
Y_A = np.cos(phi) * r
X_B = l + np.sin(phi) * r
Y_B = l + np.cos(phi) * r

fig = plt.figure(figsize=[13, 9])
ax = fig.add_subplot(1, 2, 1)
ax.axis('equal')
ax.set(xlim=[-20, 20], ylim=[-15, 15])

circle1 = plt.Circle((0, 0), r, color='black', fill=False)
ax.add_artist(circle1)

Point_A = ax.plot(X_A[0], Y_A[0], color='black')[0]
Point_B = ax.plot(X_B[0], Y_B[0], marker='o', markersize=12, color='black')[0]
Line_AB = ax.plot([X_A[0], X_B[0]], [Y_A[0], Y_B[0]], color='black', linewidth=5)[0]
#Line_OC = ax.plot([X_O, X_C[0]], [Y_O, Y_C[0]], color='black', )[0]
triangle, = ax.plot([math.sqrt(e / 2) - 0.3, math.sqrt(e / 2), math.sqrt(e / 2) + 0.3],
                    [-math.sqrt(e / 2) - 0.5, -math.sqrt(e / 2), -math.sqrt(e / 2) - 0.5], color='black')
line_tr = ax.plot([math.sqrt(e / 2) - 0.3, math.sqrt(e / 2) + 0.3], [-math.sqrt(e / 2) - 0.5, -math.sqrt(e / 2) - 0.5],
                  color='black')[0]
line_OC = ax.plot([X_C,X_O], [Y_C, Y_O], color = 'black')

# spiral spring
Nv = 1.1
R1 = 0.2
R2 = 4
thetta = np.linspace(0, Nv * 6.28 - tau[0], 100)
X_SpiralSpr = -(R1 * thetta * (R2 - R1) / thetta[-1]) * np.sin(thetta)
Y_SpiralSpr = (R1 * thetta * (R2 - R1) / thetta[-1]) * np.cos(thetta)
Drawed_Spiral_Spring = ax.plot(X_SpiralSpr + X_O, Y_SpiralSpr + Y_O, color='black')[0]

# plots
VXB = np.diff(X_B)
VYB = np.diff(Y_B)
WXB = np.diff(VXB)
WYB = np.diff(VYB)

ax2 = fig.add_subplot(4, 2, 2)
ax2.plot(VXB)
plt.title('Vx of ball')
plt.xlabel('t values')
plt.ylabel('Vx values')

ax3 = fig.add_subplot(4, 2, 4)
ax3.plot(VYB)
plt.title('Vy of ball')
plt.xlabel('t values')
plt.ylabel('Vy values')

ax4 = fig.add_subplot(4, 2, 6)
ax4.plot(WXB)
plt.title('Wx of ball')
plt.xlabel('t values')
plt.ylabel('Wy values')

ax5 = fig.add_subplot(4, 2, 8)
ax5.plot(WYB)
plt.title('Wy of ball')
plt.xlabel('t values')
plt.ylabel('Wx values')

plt.subplots_adjust(wspace=0.3, hspace=0.7)


def Dordge(i):
    Point_A.set_data(X_A[i], Y_A[i])
    Point_B.set_data(X_B[i], Y_B[i])
    Line_AB.set_data([X_A[i], X_B[i]], [Y_A[i], Y_B[i]])

    thetta = np.linspace(0, Nv * 6.28 - phi[i], 100)
    X_SpiralSpr = -(R1 * thetta * (R2 - R1) / thetta[-1]) * np.sin(thetta)
    Y_SpiralSpr = (R1 * thetta * (R2 - R1) / thetta[-1]) * np.cos(thetta)
    Drawed_Spiral_Spring.set_data(X_SpiralSpr, Y_SpiralSpr)
    #Line_OC.set_data([X_O, X_C[i]], [Y_O, Y_C[i]])
    return [Point_A, Point_B, Line_AB, Drawed_Spiral_Spring]


anim = FuncAnimation(fig, Dordge, frames=1000, interval=10)

plt.show()