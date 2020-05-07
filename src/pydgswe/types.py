from dataclasses import dataclass, field
from typing import List
from math import sqrt
import numpy as np

def initialise(geometry, f):
    p_interfaces = [f(x) for x in geometry.interfaces()]
    return [Plane(const=0.5*(p_w+p_e), slope=0.5/sqrt(3.0)*(p_e-p_w))
            for p_w, p_e in zip(p_interfaces, p_interfaces[1:])]

def piecewise(planes):
    ps_pos = [p.pos_limit() for p in planes]
    ps_neg = [p.neg_limit() for p in planes]
    return [val for pair in zip(ps_pos, ps_neg) for val in pair]

class Geometry:
    def __init__(self, elements, extent):
        self.elements = elements
        self.extent = extent
        self.dx = (extent[1] - extent[0]) / elements

    def interfaces(self):
        return np.linspace(self.extent[0], self.extent[1], self.elements+1)

@dataclass
class Plane:
    const: float = 0.0
    slope: float = 0.0

    @staticmethod
    def zero():
        return Plane()

    def neg_limit(self):
        return self.const + sqrt(3)*self.slope

    def pos_limit(self):
        return self.const - sqrt(3)*self.slope

    def gauss_west(self):
        return self.const - self.slope

    def gauss_east(self):
        return self.const + self.slope

    def __add__(self, other):
        return Plane(self.const + other.const, self.slope + other.slope)

    def __sub__(self, other):
        return Plane(self.const - other.const, self.slope - other.slope)

    def __mul__(self, scalar):
        return Plane(scalar * self.const, scalar * self.slope)

    def __rmul__(self, scalar):
        return Plane(scalar * self.const, scalar * self.slope)

    def __truediv__(self, scalar):
        return Plane(self.const / scalar, self.slope / scalar)

@dataclass
class FlowVector:
    h: float = 0.0
    q: float = 0.0

    @staticmethod
    def zero():
        return FlowVector()

    def __neg__(self):
        return FlowVector(-self.h, -self.q)

    def __add__(self, other):
        return FlowVector(self.h + other.h, self.q + other.q)

    def __sub__(self, other):
        return FlowVector(self.h - other.h, self.q - other.q)

    def __rmul__(self, scalar):
        return FlowVector(scalar * self.h, scalar * self.q)

    def __truediv__(self, scalar):
        return FlowVector(self.h / scalar, self.q / scalar)

@dataclass
class FlowCoeffs:
    h: Plane = field(default_factory=Plane)
    q: Plane = field(default_factory=Plane)

    @staticmethod
    def zero():
        return FlowCoeffs()

    def const(self):
        return FlowVector(self.h.const, self.q.const)

    def set_const(self, flow_vector):
        self.h.const = flow_vector.h
        self.q.const = flow_vector.q

    def set_slope(self, flow_vector):
        self.h.slope = flow_vector.h
        self.q.slope = flow_vector.q

    def neg_limit(self):
        return FlowVector(self.h.neg_limit(), self.q.neg_limit())

    def pos_limit(self):
        return FlowVector(self.h.pos_limit(), self.q.pos_limit())

    def gauss_west(self):
        return FlowVector(self.h.gauss_west(), self.q.gauss_west())

    def gauss_east(self):
        return FlowVector(self.h.gauss_east(), self.q.gauss_east())

    def __add__(self, other):
        return FlowCoeffs(self.h + other.h, self.q + other.q)

    def __rmul__(self, scalar):
        return FlowCoeffs(scalar * self.h, scalar * self.q)

@dataclass
class State:
    Us: List[FlowCoeffs]

    @staticmethod
    def zeros(geometry):
        return State([FlowCoeffs() for _ in range(geometry.elements)])

    @staticmethod
    def initialise(geometry, h, q = None):
        if q is None:
            q = lambda x: 0.0

        hs = initialise(geometry, h)
        qs = initialise(geometry, q)

        return State([FlowCoeffs(h, q) for h, q in zip(hs, qs)])

    def piecewise(self):
        return piecewise(self)

    def total_mass(self):
        return sum([U.h.const for U in self])

    def __add__(self, other):
        return State([U_a + U_b for U_a, U_b in zip(self.Us, other.Us)])

    def __rmul__(self, scalar):
        return State([scalar * U for U in self.Us])

    def __getitem__(self, key):
        return self.Us[key]

    def __len__(self):
        return len(self.Us)

class DEM:
    def __init__(self, zs):
        self.zs = zs
        self.zstars = [zs[0].pos_limit()]
        self.zstars += [max(z_l.neg_limit(), z_r.pos_limit())
                for z_l, z_r in zip(zs, zs[1:])]
        self.zstars += [zs[-1].neg_limit()]

    @staticmethod
    def zeros(geometry):
        return DEM([Plane() for _ in range(geometry.elements)])

    @staticmethod
    def initialise(geometry, f):
        return DEM(initialise(geometry, f))

    def piecewise(self):
        return piecewise(self.zs)
