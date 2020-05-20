from math import sqrt

def signe(x):
    return 0.0 if x==0.0 else x/abs(x)

def minmod(a1, a2, a3):
    if signe(a1) == signe(a2) and signe(a1) == signe(a3):
        return signe(a1)*min(abs(a1), abs(a2), abs(a3))
    else:
        return 0.0

def slope_limiting(u0_loc, u0_bwd, u0_fwd, u1_loc, u1_bwd, u1_fwd,
        minh, dx, emsmall, Krivo_threshold, hmin_limiter):
    if minh > hmin_limiter:
        hh = abs(dx/2.0)
        u_neg_loc = u0_loc - sqrt(3.0)*u1_loc
        u_pos_loc = u0_loc + sqrt(3.0)*u1_loc

        u_pos_bwd = u0_bwd + sqrt(3.0)*u1_bwd
        u_neg_fwd = u0_fwd - sqrt(3.0)*u1_fwd

        Ineg = abs(u_neg_loc - u_pos_bwd)
        Ipos = abs(u_neg_fwd - u_pos_loc)

        norm_h = max(abs(u0_loc - u1_loc), abs(u0_loc+u1_loc))

        DS_neg = 0.0
        DS_pos = 0.0

        if abs(norm_h) > emsmall:
            DS_neg = Ineg/(hh*norm_h)
            DS_pos = Ipos/(hh*norm_h)

        alpha = 1.0
        if max(DS_neg, DS_pos) > Krivo_threshold:
            return minmod(u1_loc, alpha*(u0_loc-u0_bwd), alpha*(u0_fwd-u0_loc))
    
    return u1_loc


