from . import FlowVector

class HLL:
    def __init__(self, physics):
        self.dry = physics.dry
        self.velocity = physics.velocity
        self.celerity = physics.celerity
        self.g = physics.g
        self.physical_flux = physics.flux

    def __call__(self, U_neg, U_pos):
        if all([self.dry(U) for U in (U_neg, U_pos)]):
            return FlowVector.zero()

        u_neg = self.velocity(U_neg)
        u_pos = self.velocity(U_pos)
        a_neg = self.celerity(U_neg)
        a_pos = self.celerity(U_pos)
        h_star = pow(0.5*(a_neg + a_pos) + 0.25*(u_neg - u_pos), 2.0) / self.g
        u_star = 0.5*(u_neg + u_pos) + a_neg - a_pos
        a_star = self.celerity(h_star)
        
        s_neg = u_pos - 2.0*a_pos \
                if self.dry(U_neg) else min(u_neg - a_neg, u_star - a_star)
        s_pos = u_neg + 2.0*a_neg \
                if self.dry(U_pos) else max(u_pos + a_pos, u_star + a_star)

        if s_neg >= 0.0:
            return self.physical_flux(U_neg)
        elif s_pos >= 0.0:
            F_neg = self.physical_flux(U_neg)
            F_pos = self.physical_flux(U_pos)

            return (s_pos*F_neg - s_neg*F_pos + s_neg*s_pos*(U_pos - U_neg)) \
                / (s_pos - s_neg)
        else:
            return self.physical_flux(U_pos)

