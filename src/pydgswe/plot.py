import matplotlib.pyplot as plt
import numpy as np

class Plot:
    def __init__(self, geometry, dem):
        self.dx = geometry.dx
        self.elements = geometry.elements
        self.xs = geometry.interfaces()
        self.xs = np.repeat(self.xs, 2)
        self.xs = self.xs[1:-1]
        self.zs = dem.piecewise()

        plt.ion()
        plt.show()

    def __call__(self, state):
        Us = state.piecewise()

        plt.clf()
        plt.subplot(311)
        plt.plot(self.xs, [z for z in self.zs])
        plt.plot(self.xs, [z+U.h for z, U in zip(self.zs, Us)])
        plt.subplot(312)
        plt.plot(self.xs, [U.h for U in Us])
        plt.subplot(313)
        plt.plot(self.xs, [U.q for U in Us])
        plt.pause(1e-4)

    def block(self):
        plt.show(block=True)