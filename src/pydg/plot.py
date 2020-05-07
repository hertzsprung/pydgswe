import matplotlib.pyplot as plt
import numpy as np

class Plot:
    def __init__(self, geometry):
        self.dx = geometry.dx
        self.elements = geometry.elements
        plt.ion()
        plt.show()

    def __call__(self, state):
        xs = np.linspace(0.0, self.dx*self.elements, self.elements+1)
        xs = np.repeat(xs, 2)
        xs = xs[1:-1]

        Us = state.piecewise()

        plt.clf()
        plt.subplot(211)
        plt.plot(xs, [U.h for U in Us])
        plt.subplot(212)
        plt.plot(xs, [U.q for U in Us])
        plt.pause(1e-4)

    def block(self):
        plt.show(block=True)
