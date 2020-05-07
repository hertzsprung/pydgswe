from dataclasses import dataclass, field
from typing import List
from math import sqrt

@dataclass
class Geometry:
    elements: int
    dx: float

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

    def __rmul__(self, scalar):
        return Plane(scalar * self.const, scalar * self.slope)

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

    def piecewise(self):
        U_poss = [U.pos_limit() for U in self]
        U_negs = [U.neg_limit() for U in self]

        return [val for pair in zip(U_poss, U_negs) for val in pair]

    def __add__(self, other):
        return State([U_a + U_b for U_a, U_b in zip(self.Us, other.Us)])

    def __rmul__(self, scalar):
        return State([scalar * U for U in self.Us])

    def __getitem__(self, key):
        return self.Us[key]

    def __len__(self):
        return len(self.Us)

@dataclass
class DEM:
    zs: List[Plane]
    zstars: List[float]

    @staticmethod
    def zeros(geometry):
        return DEM([Plane() for _ in range(geometry.elements)],
                [0.0 for _ in range(geometry.elements+1)])

