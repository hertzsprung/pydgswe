import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

def piecewise(p0, p1):
    ps_pos = [p0_ - sqrt(3.0)*p1_ for p0_, p1_ in zip(p0, p1)]
    ps_neg = [p0_ + sqrt(3.0)*p1_ for p0_, p1_ in zip(p0, p1)]
    return [val for pair in zip(ps_pos, ps_neg) for val in pair]

def plot(xx, x_interface, z0, z1, h0, h1, q0, q1, tolh):
    x_piecewise = [x_interface[0]]
    x_piecewise += list(np.repeat(x_interface[1:-1], 2))
    x_piecewise += [x_interface[-1]]

    plt.clf()
    plt.subplot(411)
    plt.gca().set_title("Elevation")
    plt.plot(x_piecewise, piecewise(z0, z1), color='C2')
    plt.plot(x_piecewise, [z+h for z, h in zip(piecewise(z0, z1), piecewise(h0, h1))], color='C0')
    plt.subplot(412)
    plt.gca().set_title("Depth")
    plt.plot(x_piecewise, piecewise(h0, h1))
    plt.subplot(413)
    plt.gca().set_title("Discharge")
    plt.plot(x_piecewise, piecewise(q0, q1))

    plt.subplot(414)
    plt.gca().set_title("Velocity")
    def velocity(h, q):
        return 0.0 if h <= tolh else q/h
    plt.plot(xx, [velocity(h, q) for h, q in zip(h0, q0)])

    plt.pause(1e-4)
