from . import DEM, Geometry, Plane, State
from math import exp, sin, cos, sqrt

class DamBreak:
    def __init__(self, elements=32):
        self.end_time = 10.0
        self.dt = 0.01
        self.geometry = Geometry(elements, extent=[0.0, 25.0])
        self.dem = DEM.zeros(self.geometry)
        self.state = State.initialise(self.geometry, self.h)

    def h(self, x):
        return 0.0 if x > 12.5 else 6.0

class LakeAtRest:
    def __init__(self, elements=32, multiplier=10.0):
        self.geometry = Geometry(elements, extent=[0.0, 50.0])
        self.end_time = 100.0
        self.dt = 0.1
        self.multiplier = multiplier
        self.dem = DEM.initialise(self.geometry, self.z)
        self.state = State.initialise(self.geometry, self.h)

    def z(self, x):
        value = 0.0

        if x >= 22.0 and x < 25.0:
            value = 0.05*x - 1.1
        elif x >= 25.0 and x <= 28.0:
            value = -0.05*x + 1.4
        elif x > 8.0 and x < 12.0:
            value = 0.2 - (0.05*(x-10.0)**2)
        elif x > 39.0 and x < 46.5:
            value =0.3

        return self.multiplier * value

    def h(self, x):
        return 2.0 - self.z(x)

class ParabolicBowlLiangMarche:
    def __init__(self, physics, elements=64, h0=10.0, a=3000.0):
        self.g = physics.g
        self.h0 = h0
        self.a = a
        self.tau = 0.0
        self.B = 8.0
        self.rho = sqrt(8.0*self.g*self.h0) / self.a

        self.end_time = 24192
        self.dt = 1.0
        self.geometry = Geometry(elements, extent=[-5000.0, 5000.0])
        self.dem = DEM.initialise(self.geometry, self.z)
        self.state = State.initialise(self.geometry, self.h)

    def z(self, x):
        return self.h0 * (x/self.a)**2

    def h(self, x):
        return max(0.0, self.eta(x) - self.z(x))

    def eta(self, x):
        t = 0.0
        s = sqrt(self.rho**2.0 - self.tau**2.0) / 2.0

        return self.h0 + (self.a**2.0) * self.B**2.0 * \
                exp(-t*self.tau) / (8.0 * self.g**2.0 * self.h0) * \
                (-s * self.tau * sin(2.0*s*t) + (self.tau**2.0 / 4.0 - s**2.0)
                        * cos(2.0*s*t)) - self.B**2.0 * \
                exp(-t*self.tau) / (4.0*self.g) - exp(-t*self.tau) / \
                (2.0*self.g) * (self.B*s*cos(s*t) + self.tau*self.B/2.0 * \
                sin(s*t)) * x

class ThinFlow:
    def __init__(self, h=1.5e-3):
        self.end_time = 120.0
        self.dt = 0.1
        self.geometry = Geometry(elements=50, extent=[0.0, 1000.0])
        self.state = State.initialise(self.geometry, lambda x: h)

        consts = [291.263908, 291.036728, 290.682968, 290.235313, 289.719070, 289.132965, 288.480148, 287.720306, 286.741554, 285.469368, 284.042336, 282.607025, 281.165619, 279.731247, 278.462189, 277.765472, 275.571877, 275.075310, 278.937813, 281.724838, 282.658119, 283.701721, 285.102501, 286.878281, 288.737343, 289.880310, 290.055466, 289.889687, 289.709999, 289.529999, 289.366722, 289.350624, 289.033592, 288.105789, 286.947197, 286.085632, 286.443604, 287.413918, 288.001411, 288.341400, 288.556717, 289.018135, 289.782349, 290.676094, 291.835312, 293.409065, 295.505478, 298.112350, 301.087814, 302.654373]
        slopes = [0.000000, -0.084529, -0.119715, -0.138739, -0.159314, -0.179074, -0.197830, -0.240865, -0.324218, -0.410279, -0.413618, -0.415059, -0.417138, -0.410997, -0.321694, -0.080556, -1.185918, 0.899225, 1.330793, 0.278297, 0.260533, 0.341991, 0.466749, 0.558498, 0.514832, 0.145060, -0.043934, -0.051779, -0.051964, -0.051959, -0.042308, 0.033014, -0.216053, -0.319615, -0.349299, -0.148126, 0.354800, 0.205411, 0.133779, 0.062513, 0.061800, 0.204600, 0.236619, 0.279385, 0.389889, 0.518718, 0.691647, 0.813432, 0.904453, 0.000000]

        self.dem = DEM([Plane(const, slope)
            for const, slope in zip(consts, slopes)])
