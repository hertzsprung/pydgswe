from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from . import *

def Baseline_Mesh(xmin, xmax, elements):
    x_int = np.linspace(xmin, xmax, elements+1)
    x_centres = [0.5*(x_w+x_e) for x_w, x_e in zip(x_int, x_int[1:])]
    return (x_centres, x_int)

def LFV(ic, _, ndir, z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp):
    z0 = z0_temp[ic]
    z1 = z1_temp[ic]
    h0 = h0_temp[ic]
    h1 = h1_temp[ic]
    q0 = q0_temp[ic]
    q1 = q1_temp[ic]

    # Eastern face
    if ndir == 2:
        et_face = (h0+z0) + sqrt(3.0)*(z1+h1)
        q_face = q0 + sqrt(3.0)*q1
        h_face = h0 + sqrt(3.0)*h1
        return (h_face, et_face, q_face)
    # Western face
    elif ndir == 4:
        et_face = (h0+z0) - sqrt(3.0)*(z1+h1)
        q_face = q0 - sqrt(3.0)*q1
        h_face = h0 - sqrt(3.0)*h1
        return (h_face, et_face, q_face)

def Wetting_Drying_1D(h_L, h_R, et_L, et_R, q_L, q_R, ndir, tolh):
    z_L = et_L - h_L
    z_R = et_R - h_R

    u_L = 0.0 if h_L <= tolh else q_L/h_L
    u_R = 0.0 if h_R <= tolh else q_R/h_R

    z_LR = max(z_L, z_R) # since z could be discontinuous at interface

    delta = 0.0
    if ndir == 2:
        delta = max(0.0, -(et_L - z_LR))
    elif ndir == 4:
        delta = max(0.0, -(et_R - z_LR))
    else:
        raise RuntimeError('invalid ndir value')

    h_L_star = max(0.0, et_L - z_LR)
    q_L_star = u_L * h_L_star

    h_R_star = max(0.0, et_R - z_LR)
    q_R_star = u_R * h_R_star

    z_LR -= delta

    return (z_LR, h_L_star, h_R_star, q_L_star, q_R_star)

def Flux_F(h, q, g, tolh):
    F = np.array([0.0, 0.0])

    if h > tolh:
        F[0] = q
        F[1] = q**2/h + g/2.0 * h**2

    return F

def DG2_1D(F_pos, F_neg, z1, h0, h1, q0, q1, dx, g, tolh):
    L0 = -(1.0/dx)*((F_pos - F_neg) + np.array([0.0, 2.0*sqrt(3.0)*g*h0*z1]))
    Flux_Q1 = Flux_F(h0-h1, q0-q1, g, tolh)
    Flux_Q2 = Flux_F(h0+h1, q0+q1, g, tolh)
    L1 = -(1.0/dx)*sqrt(3.0) * (F_pos + F_neg - (Flux_Q1 + Flux_Q2) +
            np.array([0.0, 2.0*g*h1*z1]))
    return (L0, L1)

def Flux_HLL(h_L, h_R, q_L, q_R, g, tolh):
    if h_L <= tolh and h_R <= tolh:
        return np.array([0.0, 0.0])

    u_L = 0.0 if h_L <= tolh else q_L/h_L
    u_R = 0.0 if h_R <= tolh else q_R/h_R

    a_L = sqrt(g*h_L)
    a_R = sqrt(g*h_R)

    h_star = ((a_L+a_R)/2.0 + (u_L-u_R)/4.0)**2/g
    u_star = (u_L+u_R)/2.0 + a_L - a_R
    a_star = sqrt(g*h_star)

    s_L = u_R-2*a_R if h_L <= tolh else min(u_L-a_L, u_star-a_star)
    s_R = u_L+2*a_L if h_R <= tolh else max(u_R+a_R, u_star+a_star)

    F_L = np.array([q_L, u_L*q_L+0.5*g*h_L**2])
    F_R = np.array([q_R, u_R*q_R+0.5*g*h_R**2])

    if s_L >= 0:
        return F_L
    elif s_L < 0 and s_R >= 0:
        return np.array([
            (s_R*F_L[0]-s_L*F_R[0]+s_L*s_R*(h_R-h_L))/(s_R-s_L),
            (s_R*F_L[1]-s_L*F_R[1]+s_L*s_R*(q_R-q_L))/(s_R-s_L)
        ])
    else:
        return F_R

def DG2_Operator(ic, dx, elements, h_L_up, q_L_up, h_R_dw, q_R_dw,
        z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp, g, tolh):
    # L and R denote the interface limits approaching from the left and right

    # Eastern face (pos)
    h_L, et_L, q_L = LFV(ic, dx, 2,
            z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp)

    h_R = 0.0
    et_R = 0.0
    q_R = 0.0
    if ic == elements: # eastern-most elemnt (ignoring ghosts)
        h_R = h_R_dw
        et_R = (et_L-h_L)+h_R
        q_R = q_R_dw
    else:
        h_R, et_R, q_R = LFV(ic+1, dx, 4, 
            z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp)

    zpos_LR, h_pos_L_star, h_pos_R_star, q_pos_L_star, q_pos_R_star = \
            Wetting_Drying_1D(h_L, h_R, et_L, et_R, q_L, q_R, 2, tolh)
    
    F_pos = Flux_HLL(h_pos_L_star, h_pos_R_star, q_pos_L_star, q_pos_R_star,
            g, tolh)
    
    # Western face (neg)
    h_R, et_R, q_R = LFV(ic, dx, 4,
            z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp)

    h_L = 0.0
    et_L = 0.0
    q_L = 0.0
    if ic == 1: # western-most element (ignoring ghosts)
        h_L = h_L_up
        et_L = (et_R-h_R)+h_L
        q_L = q_L_up
    else:
        h_L, et_L, q_L = LFV(ic-1, dx, 2, 
            z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp)

    zneg_LR, h_neg_L_star, h_neg_R_star, q_neg_L_star, q_neg_R_star = \
            Wetting_Drying_1D(h_L, h_R, et_L, et_R, q_L, q_R, 4, tolh)

    F_neg = Flux_HLL(h_neg_L_star, h_neg_R_star, q_neg_L_star, q_neg_R_star,
            g, tolh)

    # positivity preserving coefficients (denoted bar)
    z1_bar = (zpos_LR - zneg_LR)*(sqrt(3.0)/6.0)
    h0_bar = (h_pos_L_star + h_neg_R_star) / 2.0
    h1_bar = (h_pos_L_star - h_neg_R_star)*(sqrt(3.0)/6.0)
    q0_bar = (q_pos_L_star + q_neg_R_star) / 2.0
    q1_bar = (q_pos_L_star - q_neg_R_star)*(sqrt(3.0)/6.0)

    return DG2_1D(F_pos, F_neg, z1_bar, h0_bar, h1_bar, q0_bar, q1_bar, dx,
            g, tolh)

def friction_implicit(h0, h1, q0, q1, dt, g, tolh, manning, emsmall):
    h_G1 = h0 + h1
    h_G2 = h0 - h1
    q_G1 = q0 - q1
    q_G2 = q0 + q1

    minh = min(h_G1, h0, h_G2)
    minq = min(q_G1, q0, q_G2)

    q0_friction = q0
    q1_friction = q1

    if minh > tolh and abs(minq) > emsmall:
        u0 = q0/h0
        Cf0 = g*manning**2/pow(h0, 1.0/3.0)
        Sf0 = -Cf0*abs(u0)*u0
        D0 = 1.0 + 2.0*dt*Cf0*abs(u0)/h0
        q0_friction = q0 + dt*Sf0/D0

        u_G1 = q_G1/h_G1
        Cf_G1 = g*manning**2/pow(h_G1, 1.0/3.0)
        Sf_G1 = -Cf_G1*abs(u_G1)*u_G1
        D_G1 = 1.0 + 2.0*dt*Cf_G1*abs(u_G1)/h_G1
        q_G1_friction = q_G1 + dt*Sf_G1/D_G1

        u_G2 = q_G2/h_G2
        Cf_G2 = g*manning**2/pow(h_G2, 1.0/3.0)
        Sf_G2 = -Cf_G2*abs(u_G2)*u_G2
        D_G2 = 1.0 + 2.0*dt*Cf_G2*abs(u_G2)/h_G2
        q_G2_friction = q_G2 + dt*Sf_G2/D_G2

        q1_friction = 1.0/2.0 * (q_G1_friction - q_G2_friction)
        
    return (q0_friction, q1_friction)

def main():
    tolh = 1e-3
    emsmall = 1e-12
    time_now = 0.0
    CFL = 0.28
    Krivo_threshold = 10.0

#    case = LakeAtRest()
#    case = BuildingOvertopping()
    case = SheetFlow()

    xmin = case.xmin
    xmax = case.xmax
    elements = 50
    bcs = case.bcs

    simulation_time = case.simulation_time
    g = 9.80665
    manning = case.manning

    dx = (xmax - xmin) / elements
    xx, x_interface = Baseline_Mesh(xmin, xmax, elements)

    z_interface = [case.Bed_Data(x) for x in x_interface] 
    h_interface = [case.Init_Conds_h(z, x)
            for z, x in zip(z_interface, x_interface)]
    q_interface = [case.Init_Conds_q(z, x)
            for z, x in zip(z_interface, x_interface)]

    z0 = [0.5*(z_w+z_e) for z_w, z_e in zip(z_interface, z_interface[1:])]
    h0 = [0.5*(h_w+h_e) for h_w, h_e in zip(h_interface, h_interface[1:])]
    q0 = [0.5*(q_w+q_e) for q_w, q_e in zip(q_interface, q_interface[1:])]
    
    z1 = [(z_e-z_w)*(sqrt(3.0)/6.0)
            for z_w, z_e in zip(z_interface, z_interface[1:])]
    h1 = [(h_e-h_w)*(sqrt(3.0)/6.0)
            for h_w, h_e in zip(h_interface, h_interface[1:])]
    q1 = [(q_e-q_w)*(sqrt(3.0)/6.0)
            for q_w, q_e in zip(q_interface, q_interface[1:])]

    h0_max = max(h0)
    hmin_limiter = h0_max*5.0/100.0
    dt = 0.0001

    plt.ion()
    plt.show()
    plot(xx, x_interface, z0, z1, h0, h1, q0, q1, tolh)

    c = 0
    while time_now < simulation_time:
        time_now += dt
        if time_now - simulation_time > 0:
            time_now -= dt
            dt = simulation_time - time_now
            time_now += dt

        print("t", time_now, "dt", dt)

        h0_with_bc, z0_with_bc, q0_with_bc, \
                h1_with_bc, z1_with_bc, q1_with_bc, \
                h_L_up, q_L_up, h_R_dw, q_R_dw \
                = bcs.Add_Ghost_BCs(h0, z0, q0, h1, z1, q1)

        # intermediate data after RK stage 1
        h0_int_with_bc = h0_with_bc.copy()
        h1_int_with_bc = h1_with_bc.copy()
        q0_int_with_bc = q0_with_bc.copy()
        q1_int_with_bc = q1_with_bc.copy()

        # new data after RK stage 2
        h0_new_with_bc = h0_with_bc.copy()
        h1_new_with_bc = h1_with_bc.copy()
        q0_new_with_bc = q0_with_bc.copy()
        q1_new_with_bc = q1_with_bc.copy()

        if manning > 0.0:
            for i in range(len(h0_with_bc)):
                q0_friction, q1_friction = friction_implicit(
                        h0_with_bc[i], h1_with_bc[i],
                        q0_with_bc[i], q1_with_bc[i], dt,
                        g, tolh, manning, emsmall)
                q0_with_bc[i] = q0_friction
                q1_with_bc[i] = q1_friction

        # RK stage 1
        z0_temp = z0_with_bc.copy()
        z1_temp = z1_with_bc.copy()
        h0_temp = h0_with_bc.copy()
        h1_temp = h1_with_bc.copy()
        q0_temp = q0_with_bc.copy()
        q1_temp = q1_with_bc.copy()

        # slope limiting
        for i in range(1, elements+1):
            et0_loc = h0_temp[i] + z0_temp[i]
            et1_loc = h1_temp[i] + z1_temp[i]
            et0_bwd = h0_temp[i-1] + z0_temp[i-1]
            et1_bwd = h1_temp[i-1] + z1_temp[i-1]
            et0_fwd = h0_temp[i+1] + z0_temp[i+1]
            et1_fwd = h1_temp[i+1] + z1_temp[i+1]
            q0_loc = q0_temp[i]
            q1_loc = q1_temp[i]
            q0_bwd = q0_temp[i-1]
            q1_bwd = q1_temp[i-1]
            q0_fwd = q0_temp[i+1]
            q1_fwd = q1_temp[i+1]

            h0_loc = h0_temp[i]
            h_G1_loc = h0_temp[i] + h1_temp[i]
            h_G2_loc = h0_temp[i] - h1_temp[i]
            min_h_loc = min(h0_loc, h_G1_loc, h_G2_loc)

            h0_bwd = h0_temp[i-1]
            h_G1_bwd = h0_temp[i-1]+h1_temp[i-1]
            h_G2_bwd = h0_temp[i-1]-h1_temp[i-1]
            min_h_bwd = min(h0_bwd, h_G1_bwd, h_G2_bwd)

            h0_fwd = h0_temp[i+1]
            h_G1_fwd = h0_temp[i+1]+h1_temp[i+1]
            h_G2_fwd = h0_temp[i+1]-h1_temp[i+1]
            min_h_fwd = min(h0_fwd, h_G1_fwd, h_G2_fwd)

            min_h = min(min_h_bwd, min_h_loc, min_h_fwd)

            et1_loc = slope_limiting(et0_loc, et0_bwd, et0_fwd, et1_loc, et1_bwd, et1_fwd, min_h, dx, emsmall, Krivo_threshold, hmin_limiter)
            q1_loc = slope_limiting(q0_loc, q0_bwd, q0_fwd, q1_loc, q1_bwd, q1_fwd, min_h, dx, emsmall, Krivo_threshold, hmin_limiter)

            h1_temp[i] = et1_loc - z1_temp[i]
            q1_temp[i] = q1_loc

            h1_with_bc[i] = h1_temp[i]
            q1_with_bc[i] = q1_temp[i]
        

        for i in range(1, elements+1):
            L0_temp, L1_temp = DG2_Operator(i, dx, elements,
                    h_L_up, q_L_up, h_R_dw, q_R_dw,
                    z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp,
                    g, tolh)

            h0_int_with_bc[i] = h0_with_bc[i] + dt*L0_temp[0]
            q0_int_with_bc[i] = q0_with_bc[i] + dt*L0_temp[1]
            h1_int_with_bc[i] = h1_with_bc[i] + dt*L1_temp[0]
            q1_int_with_bc[i] = q1_with_bc[i] + dt*L1_temp[1]

            if h0_int_with_bc[i] <= tolh:
                q0_int_with_bc[i] = 0.0
                q1_int_with_bc[i] = 0.0

        # RK stage 2
        h0_int_with_bc, _, q0_int_with_bc, \
                h1_int_with_bc, _, q1_int_with_bc, \
                h_L_up, q_L_up, h_R_dw, q_R_dw \
                = bcs.Add_Ghost_BCs(
                        h0_int_with_bc[1:-1],
                        z0,
                        q0_int_with_bc[1:-1],
                        h1_int_with_bc[1:-1],
                        z1,
                        q1_int_with_bc[1:-1])

        h0_temp = h0_int_with_bc
        h1_temp = h1_int_with_bc
        q0_temp = q0_int_with_bc
        q1_temp = q1_int_with_bc

        # slope limiting
        for i in range(1, elements+1):
            et0_loc = h0_temp[i] + z0_temp[i]
            et1_loc = h1_temp[i] + z1_temp[i]
            et0_bwd = h0_temp[i-1] + z0_temp[i-1]
            et1_bwd = h1_temp[i-1] + z1_temp[i-1]
            et0_fwd = h0_temp[i+1] + z0_temp[i+1]
            et1_fwd = h1_temp[i+1] + z1_temp[i+1]
            q0_loc = q0_temp[i]
            q1_loc = q1_temp[i]
            q0_bwd = q0_temp[i-1]
            q1_bwd = q1_temp[i-1]
            q0_fwd = q0_temp[i+1]
            q1_fwd = q1_temp[i+1]

            h0_loc = h0_temp[i]
            h_G1_loc = h0_temp[i] + h1_temp[i]
            h_G2_loc = h0_temp[i] - h1_temp[i]
            min_h_loc = min(h0_loc, h_G1_loc, h_G2_loc)

            h0_bwd = h0_temp[i-1]
            h_G1_bwd = h0_temp[i-1]+h1_temp[i-1]
            h_G2_bwd = h0_temp[i-1]-h1_temp[i-1]
            min_h_bwd = min(h0_bwd, h_G1_bwd, h_G2_bwd)

            h0_fwd = h0_temp[i+1]
            h_G1_fwd = h0_temp[i+1]+h1_temp[i+1]
            h_G2_fwd = h0_temp[i+1]-h1_temp[i+1]
            min_h_fwd = min(h0_fwd, h_G1_fwd, h_G2_fwd)

            min_h = min(min_h_bwd, min_h_loc, min_h_fwd)

            et1_loc = slope_limiting(et0_loc, et0_bwd, et0_fwd, et1_loc, et1_bwd, et1_fwd, min_h, dx, emsmall, Krivo_threshold, hmin_limiter)
            q1_loc = slope_limiting(q0_loc, q0_bwd, q0_fwd, q1_loc, q1_bwd, q1_fwd, min_h, dx, emsmall, Krivo_threshold, hmin_limiter)

            h1_temp[i] = et1_loc - z1_temp[i]
            q1_temp[i] = q1_loc

            h1_int_with_bc[i] = h1_temp[i]
            q1_int_with_bc[i] = q1_temp[i]

        for i in range(1, elements+1):
            L0_temp, L1_temp = DG2_Operator(i, dx, elements,
                    h_L_up, q_L_up, h_R_dw, q_R_dw,
                    z0_temp, z1_temp, h0_temp, h1_temp, q0_temp, q1_temp,
                    g, tolh)

            h0_new_with_bc[i] = 0.5*(h0_with_bc[i] + h0_int_with_bc[i]
                    + dt*L0_temp[0])
            q0_new_with_bc[i] = 0.5*(q0_with_bc[i] + q0_int_with_bc[i]
                    + dt*L0_temp[1])
            h1_new_with_bc[i] = 0.5*(h1_with_bc[i] + h1_int_with_bc[i]
                    + dt*L1_temp[0])
            q1_new_with_bc[i] = 0.5*(q1_with_bc[i] + q1_int_with_bc[i]
                    + dt*L1_temp[1])

            if h0_new_with_bc[i] <= tolh:
                q0_new_with_bc[i] = 0.0
                q1_new_with_bc[i] = 0.0

        # calculate next timestep
        dt = 1e9
        for i in range(elements):
            h0[i] = h0_new_with_bc[i+1]
            q0[i] = q0_new_with_bc[i+1]
            h1[i] = h1_new_with_bc[i+1]
            q1[i] = q1_new_with_bc[i+1]

            h_G1 = h0[i] - h1[i]
            h_G2 = h0[i] + h1[i]
            q_G1 = q0[i] - q1[i]
            q_G2 = q0[i] + q1[i]

            if h0[i] <= tolh:
                q0[i] = 0.0
                q1[i] = 0.0
                continue
            else:
                if h0[i] > tolh:
                    u0 = q0[i]/h0[i]
                    dt0 = CFL*dx/(abs(u0)+sqrt(g*h0[i]))
                    dt = min(dt, dt0)
                #if h_G1 > tolh:
                #    u_G1 = q_G1/h_G1
                #    dt_G1 = CFL*dx/(abs(u_G1)+sqrt(g*h_G1))
                #    dt = min(dt, dt_G1)

                #if h_G2 > tolh:
                #    u_G2 = q_G2/h_G2
                #    dt_G2 = CFL*dx/(abs(u_G2)+sqrt(g*h_G2))
                #    dt = min(dt, dt_G2)

        h0_max = max(h0)
        hmin_limiter = h0_max*5.0/100.0

        c += 1
        if c % 20 == 0:
            plot(xx, x_interface, z0, z1, h0, h1, q0, q1, tolh)

    plt.show(block=True)
