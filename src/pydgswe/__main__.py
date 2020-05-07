from dataclasses import dataclass
from math import sqrt
from . import DEM, FlowVector, FlowCoeffs, Geometry, HLL, State, Plane
from . import ParabolicBowlLiangMarche
from . import Plot

class RungeKutta2:
    def __init__(self, spatial_op, after_stage_op = None):
        if after_stage_op is None:
            after_stage_op = lambda state: state
            
        self.spatial_op = spatial_op
        self.after_stage_op = after_stage_op

    def __call__(self, state, dt):
        state_n = state

        state_int = state_n + dt * self.spatial_op(state)
        state_int = self.after_stage_op(state_int)

        state_new = 0.5*(state_n + state_int + dt * self.spatial_op(state_int))
        state_new = self.after_stage_op(state_new)

        return state_new

@dataclass
class Physics:
    g: float = 9.80665
    depth_threshold: float = 1e-3

    def dry(self, value):
        if isinstance(value, float):
            return value <= self.depth_threshold
        else:
            return value.h <= self.depth_threshold

    def velocity(self, U):
        return 0.0 if self.dry(U) else U.q / U.h

    def celerity(self, value):
        if isinstance(value, float):
            return sqrt(self.g * value)
        else:
            return sqrt(self.g * value.h)

    def flux_from_discharge(self, flow_vector):
        U = flow_vector
        if self.dry(U):
            return FlowVector()
        else:
            return FlowVector(U.q, U.q*U.q/U.h + 0.5*self.g*U.h*U.h)

    def flux_from_velocity(self, flow_vector):
        U = flow_vector
        return FlowVector(U.q, self.velocity(U) * U.q + 0.5*self.g*U.h*U.h)

    def star_limit(self, flow_vector, z, zstar):
        U = flow_vector
        hstar = max(0.0, z + U.h - zstar)
        qstar = 0.0 if self.dry(hstar) else hstar * self.velocity(U)
        return FlowVector(hstar, qstar)

    def star_coeffs(self, U, z, zstar_w, zstar_e):
        Ustar_pos = self.star_limit(U.pos_limit(), z.pos_limit(), zstar_w)
        Ustar_neg = self.star_limit(U.neg_limit(), z.neg_limit(), zstar_e)

        coeffs = FlowCoeffs()
        coeffs.set_const(0.5*(Ustar_pos + Ustar_neg))
        coeffs.set_slope((Ustar_neg - Ustar_pos) / (2.0*sqrt(3.0)))
        return coeffs

class DG2SpatialOperator:
    def __init__(self, riemann_solver, geometry, dem, physics):
        self.g = physics.g
        self.riemann_solver = riemann_solver
        self.dx = geometry.dx
        self.zs = dem.zs
        self.zstars = dem.zstars
        self.physical_flux = physics.flux_from_discharge
        self.star_limit = physics.star_limit
        self.star_coeffs = physics.star_coeffs

    def __call__(self, state):
        def solve(U_l, U_r, z_l, z_r, zstar):
            return self.riemann_solver(
                self.star_limit(U_l.neg_limit(), z_l.neg_limit(), zstar),
                self.star_limit(U_r.pos_limit(), z_r.pos_limit(), zstar))

        Fs = [FlowVector()] # FIXME: western BC
        Fs += [solve(U_l, U_r, z_l, z_r, zstar)
                for U_l, U_r, z_l, z_r, zstar
                in zip(state, state[1:], self.zs, self.zs[1:], self.zstars[1:])]
        Fs += [FlowVector()] # FIXME: eastern BC

        return State([self.L(U, z, zstar_w, zstar_e, F_w, F_e)
            for U, z, zstar_w, zstar_e, F_w, F_e
            in zip(state, self.zs, self.zstars, self.zstars[1:], Fs, Fs[1:])])

    def L(self, U, z, zstar_w, zstar_e, F_w, F_e):
        Ustar = self.star_coeffs(U, z, zstar_w, zstar_e)
        coeffs = FlowCoeffs()
        coeffs.set_const(-(F_e - F_w) / self.dx)
        coeffs.set_slope(-sqrt(3.0) / self.dx * (F_w + F_e
                - self.gauss_quadrature(Ustar, self.physical_flux)))
        coeffs += self.bed_slope_source(U, Ustar, z, zstar_w, zstar_e)
        return coeffs

    def gauss_quadrature(self, U, f):
        return f(U.gauss_west()) + f(U.gauss_east())

    def bed_slope_source(self, U, Ustar, z, zstar_w, zstar_e):
        zdagger_w = self.zdagger((U.h + z).pos_limit(), zstar_w)
        zdagger_e = self.zdagger((U.h + z).neg_limit(), zstar_e)
        zdagger_slope = (zdagger_e - zdagger_w) / (2.0*sqrt(3.0))

        Sb = FlowCoeffs()
        Sb.q = -2.0*sqrt(3.0) * self.g * Ustar.h * zdagger_slope / self.dx
        return Sb

    def zdagger(self, eta, zstar):
        return zstar - max(0.0, -(eta - zstar))

class ZeroDryDischarge:
    def __init__(self, physics):
        self.physics = physics

    def __call__(self, state):
        for U in state:
            if self.physics.dry(U.const()):
                U.q = Plane.zero()

        return state

def main():
    t = 0.0
    dt = 1.0
    physics = Physics()
    riemann_solver = HLL(physics)
    case = ParabolicBowlLiangMarche(physics)
    L = DG2SpatialOperator(riemann_solver, case.geometry, case.dem, physics)
    zero_dry_discharge = ZeroDryDischarge(physics)
    rk2 = RungeKutta2(L, zero_dry_discharge)

    state = case.state
    plot = Plot(case.geometry, case.dem)

    c = 0
    while t < case.end_time:
        state = rk2(state, dt)
        t += dt
        c += 1
        print("t=", t, "mass=", state.total_mass())
        if c % 50 == 0:
            plot(state)

    plot.block()
