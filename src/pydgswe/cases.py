from . import DEM, Geometry, State
from math import exp, sin, cos, sqrt

class ParabolicBowlLiangMarche:
    def __init__(self, physics, elements=64, h0=10.0, a=3000.0):
        self.g = physics.g
        self.h0 = h0
        self.a = a
        self.tau = 0.0
        self.B = 8.0
        self.rho = sqrt(8.0*self.g*self.h0) / self.a

        self.end_time = 24192
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

