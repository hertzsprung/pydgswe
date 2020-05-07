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

        h_poss = [U.h.positive_limit() for U in state]
        h_negs = [U.h.negative_limit() for U in state]

        hs = [val for pair in zip(h_poss, h_negs) for val in pair]

        plt.clf()
        plt.plot(xs, hs)
        plt.pause(1e-4)
