from dataclasses import dataclass

@dataclass
class BoundaryConditions:
    reflect_up: float = 1.0
    reflect_dw: float = 1.0
    h_imposed_up: float = 0.0
    q_imposed_up: float = 0.0
    h_imposed_dw: float = 0.0
    q_imposed_dw: float = 0.0

    def Add_Ghost_BCs(self, h0, z0, q0, h1, z1, q1):
        z0_up = z0[0]
        h0_up = h0[0]
        q0_up = self.reflect_up*q0[0]

        if self.h_imposed_up > 0.0:
            h0_up = self.h_imposed_up

        if self.q_imposed_up > 0.0:
            q0_up = self.q_imposed_up

        h1_up = 0.0
        z1_up = 0.0
        q1_up = 0.0

        z0_dw = z0[-1]
        h0_dw = h0[-1]
        q0_dw = self.reflect_dw*q0[-1]

        if self.h_imposed_dw > 0.0:
            h0_dw = self.h_imposed_dw

        if self.q_imposed_dw > 0.0:
            q0_dw = self.q_imposed_dw

        h1_dw = 0.0
        z1_dw = 0.0
        q1_dw = 0.0

        z0_with_bc = [z0_up] + z0 + [z0_dw]
        h0_with_bc = [h0_up] + h0 + [h0_dw]
        q0_with_bc = [q0_up] + q0 + [q0_dw]
        z1_with_bc = [z1_up] + z1 + [z1_dw]
        h1_with_bc = [h1_up] + h1 + [h1_dw]
        q1_with_bc = [q1_up] + q1 + [q1_dw]

        return (h0_with_bc, z0_with_bc, q0_with_bc,
                h1_with_bc, z1_with_bc, q1_with_bc,
                h0_up, q0_up, h0_dw, q0_dw)
