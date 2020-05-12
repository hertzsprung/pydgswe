from dataclasses import dataclass
from math import sqrt
from . import *

class RungeKutta2:
    def __init__(self, spatial_op, after_stage_op = None):
        if after_stage_op is None:
            after_stage_op = lambda state: state
            
        self.spatial_op = spatial_op
        self.after_stage_op = after_stage_op

    def __call__(self, state, dt):
        L = self.spatial_op
        state_n = state

        state_int = state_n + dt*L(state)
        state_int = self.after_stage_op(state_int)

        state_new = 0.5*(state_n + state_int + dt*L(state_int))
        state_new = self.after_stage_op(state_new)

        return state_new

@dataclass
class Physics:
    g: float = 9.80665
    depth_threshold: float = 1e-3
    speed_threshold: float = 1e-12
    cfl : float = 0.25
    max_dt : float = 60.0

    def dry(self, value):
        if isinstance(value, float):
            return value <= self.depth_threshold
        else:
            return value.h <= self.depth_threshold

    def motionless(self, U):
        return abs(self.velocity(U)) <= self.speed_threshold 

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

    def Ustar_at_limit(self, flow_vector, z, zstar):
        U = flow_vector
        hstar = max(0.0, z + U.h - zstar)
        qstar = 0.0 if self.dry(hstar) else hstar * self.velocity(U)
        return FlowVector(hstar, qstar)

    def Ustar_coeffs(self, U, z, zstar_w, zstar_e):
        Ustar_w = self.Ustar_at_limit(U.pos_limit(), z.pos_limit(), zstar_w)
        Ustar_e = self.Ustar_at_limit(U.neg_limit(), z.neg_limit(), zstar_e)

        return FlowCoeffs.reconstruct_from_limits(Ustar_w, Ustar_e)

    def dt(self, U, dx):
        U = U.const()

        if self.dry(U):
            return self.max_dt
        else:
            return min(self.max_dt,
                    self.cfl*dx / (abs(self.velocity(U)) + sqrt(self.g*U.h)))

class DG2SpatialOperator:
    def __init__(self, riemann_solver, geometry, dem, physics,
            boundary_condition_l, boundary_condition_r):
        self.g = physics.g
        self.riemann_solver = riemann_solver
        self.dx = geometry.dx
        self.zs = dem.zs
        self.zstars = dem.zstars
        self.physical_flux = physics.flux_from_discharge
        self.Ustar_at_limit = physics.Ustar_at_limit
        self.Ustar_coeffs = physics.Ustar_coeffs
        self.bc_l = boundary_condition_l
        self.bc_r = boundary_condition_r

    def __call__(self, state):
        def solve(U_l, U_r, z_l, z_r, zstar):
            return self.riemann_solver(
                self.Ustar_at_limit(U_l.neg_limit(), z_l.neg_limit(), zstar),
                self.Ustar_at_limit(U_r.pos_limit(), z_r.pos_limit(), zstar))

        Fs = [self.boundary_left(state)]
        Fs += [solve(U_l, U_r, z_l, z_r, zstar)
                for U_l, U_r, z_l, z_r, zstar
                in zip(state, state[1:], self.zs, self.zs[1:], self.zstars[1:])]
        Fs += [self.boundary_right(state)]

        return State([self.L(U, z, zstar_w, zstar_e, F_w, F_e)
            for U, z, zstar_w, zstar_e, F_w, F_e
            in zip(state, self.zs, self.zstars, self.zstars[1:], Fs, Fs[1:])])

    def L(self, U, z, zstar_w, zstar_e, F_w, F_e):
        Ustar = self.Ustar_coeffs(U, z, zstar_w, zstar_e)

        coeffs = FlowCoeffs()
        coeffs.set_const(-(F_e - F_w) / self.dx)
        coeffs.set_slope(-sqrt(3.0) / self.dx * (F_w + F_e
                - self.gauss_quadrature(Ustar, self.physical_flux)))
        coeffs += self.bed_slope_source(U, Ustar, z, zstar_w, zstar_e)

#        if U.h.const > 1e-3:
#            print("DEM", z, "zstar_w", zstar_w, "zstar_e",zstar_e)
#            print("h", U.h)
#            print("F_w", F_w, "F_e", F_e)
#            print("Ustar", Ustar)
#            print("flux_gradient", FlowCoeffs(-(F_e - F_w) / self.dx, -sqrt(3.0) / self.dx * (F_w + F_e
#                - self.gauss_quadrature(Ustar, self.physical_flux))))
#            print("bed_slope_source", self.bed_slope_source(U, Ustar, z, zstar_w, zstar_e))
        return coeffs

    def gauss_quadrature(self, U, f):
        return f(U.gauss_west()) + f(U.gauss_east())

    def bed_slope_source(self, U, Ustar, z, zstar_w, zstar_e):
        zdagger_w = self.zdagger((U.h + z).pos_limit(), zstar_w)
        zdagger_e = self.zdagger((U.h + z).neg_limit(), zstar_e)
        zdagger_slope = slope(zdagger_w, zdagger_e)

        coeffs = FlowCoeffs.zero()
        coeffs.q = -2.0*sqrt(3.0) * self.g * Ustar.h * zdagger_slope / self.dx
        return coeffs

    def zdagger(self, eta, zstar):
        return zstar - max(0.0, -(eta - zstar))

    def boundary_left(self, state):
        Ustar = self.Ustar_coeffs(state[0], self.zs[0],
                self.zstars[0], self.zstars[1])
        Ustar_at_limit = self.Ustar_at_limit(
                state[0].pos_limit(), self.zs[0].pos_limit(), self.zstars[0])
        return self.riemann_solver(*self.bc_l(Ustar, Ustar_at_limit))

    def boundary_right(self, state):
        Ustar = self.Ustar_coeffs(state[-1], self.zs[-1],
                self.zstars[-2], self.zstars[-1])
        Ustar_limit = self.Ustar_at_limit(
                state[-1].neg_limit(), self.zs[-1].neg_limit(), self.zstars[-1])
        return self.riemann_solver(*self.bc_r(Ustar, Ustar_limit)[::-1])

class ZeroDryDischarge:
    def __init__(self, physics):
        self.dry = physics.dry

    def __call__(self, state):
        for U in state:
            if self.dry(U.const()):
                U.q = Plane.zero()

        return state

class ZeroDryDischargeFromLimits:
    def __init__(self, physics):
        self.dry = physics.dry

    def __call__(self, state):
        for U in state:
            q_w = 0.0 if self.dry(U.pos_limit()) else U.q.pos_limit()
            q_e = 0.0 if self.dry(U.neg_limit()) else U.q.neg_limit()

            # use conditional to avoid more rounding off... but is it necessary?
            if self.dry(U.pos_limit()) or self.dry(U.neg_limit()):
                U.q = Plane.reconstruct_from_limits(q_w, q_e)

        return state

class SplittingImplicitFriction:
    def __init__(self, physics, manning):
        self.dry = physics.dry
        self.motionless = physics.motionless
        self.velocity = physics.velocity
        self.g = physics.g
        self.n = manning

    def __call__(self, state, dt):
        def friction(U):
            assert U.h >= 0.0

            u = self.velocity(U)
            Cf = self.g*self.n*self.n / pow(U.h, 1.0/3.0)
            Sf = -Cf*abs(u)*u
            D = 1.0 + 2.0*dt*Cf * abs(u) / U.h

            assert U.q * (U.q - dt*Sf/D) >= 0.0, 'Flow reversed by friction update'

            return U.q + dt*Sf/D
            
        for U in state:
            if self.dry(U.const()) or self.motionless(U.const()):
                continue
            
            q_gauss_west = friction(U.gauss_west())
            q_gauss_east = friction(U.gauss_east())

            U.q.const = friction(U.const())
            U.q.slope = 0.5*(q_gauss_east - q_gauss_west)

        return state

def main():
    t = 0.0
    physics = Physics()
    riemann_solver = HLL(physics)
    #case = ParabolicBowlLiangMarche(physics)
    #case = DamBreak()
    #case = LakeAtRest()
    case = ThinFlow()
    #case = ThinFlowSixElements()
    L = DG2SpatialOperator(riemann_solver, case.geometry, case.dem, physics,
            TransmissiveBoundary(), TransmissiveBoundary())
    zero_dry_discharge = ZeroDryDischarge(physics)
    rk2 = RungeKutta2(L, zero_dry_discharge)
    friction = lambda state, dt: state
    #friction = SplittingImplicitFriction(physics, case.manning)

    dt = case.dt
    state = case.state
    plot = Plot(case.geometry, case.dem, physics)
    plot.dts += [(t, dt)]
    plot(state)

    c = 0
    while t < case.end_time:
        dt = state.min_dt(physics, case.geometry)
        state = friction(state, dt)
        state = rk2(state, dt)
        t += dt
        c += 1
        print("t=", t,
                "mass=", state.total_mass(),
                "wet_cells=", state.total_wet(physics),
                "dry_cells=", state.total_dry(physics),
                "min_dt=", state.min_dt(physics, case.geometry))

        plot.dts += [(t, dt)]
        if c % 20 == 0:
            plot(state)

    print("total timesteps=", len(plot.dts))
    plot.block()
