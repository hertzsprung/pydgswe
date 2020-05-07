from dataclasses import dataclass
from math import sqrt
from . import Plane, FlowVector, FlowCoeffs, Geometry, HLL, State
from . import Plot

def rk2(state, op, dt):
    state_n = state
    state_int = state_n + dt * op(state)
    return 0.5 * (state_n + state_int + dt * op(state_int))

@dataclass
class Physics:
    g: float = 9.80665
    depth_threshold: float = 1e-3

    def dry(self, U):
        return U.h <= self.depth_threshold

    def velocity(self, U):
        return 0.0 if self.dry(U) else U.q / U.h

    def celerity(self, value):
        if isinstance(value, float):
            return sqrt(self.g * value)
        else:
            return sqrt(self.g * value.h)

    def flux(self, flow_vector):
        U = flow_vector
        return FlowVector(U.q, self.velocity(U) * U.q + 0.5*self.g*U.h*U.h)

class DG2SpatialOperator:
    def __init__(self, riemann_solver, geometry, physics):
        self.riemann_solver = riemann_solver
        self.dx = geometry.dx
        self.physical_flux = physics.flux

    def __call__(self, state):
        Fs = [FlowVector()] # FIXME: western BC
        Fs += [self.riemann_solver(
                    U_left.negative_limit(), U_right.positive_limit())
                for U_left, U_right in zip(state, state[1:])]
        Fs += [FlowVector()] # FIXME: eastern BC

        return State([self.L(U, F_w, F_e)
            for U, F_w, F_e in zip(state, Fs, Fs[1:])])

    def L(self, U, F_w, F_e):
        coeffs = FlowCoeffs()
        coeffs.set_const(-(F_e - F_w) / self.dx)
        coeffs.set_slope(-sqrt(3.0)/self.dx * (F_w + F_e
                - self.gauss_quadrature(U, self.physical_flux)))
        return coeffs

    def gauss_quadrature(self, U, f):
        return f(U.gauss_west()) + f(U.gauss_east())

def main():
    end_time = 20.0
    t = 0.0
    dt = 0.1
    ddt = rk2
    physics = Physics()
    riemann_solver = HLL(physics)
    geometry = Geometry(elements=40, dx=20.0)
    L = DG2SpatialOperator(riemann_solver, geometry, physics)

    state = State.zeros(geometry)
    plot = Plot(geometry)

    for i, U in enumerate(state):
        U.h.const = 5.0 if i > 10 and i < 20 else 0

    while t < end_time:
        state = ddt(state, L, dt)
        t += dt
        print(t)
        plot(state)
