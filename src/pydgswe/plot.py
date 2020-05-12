import matplotlib.pyplot as plt
import numpy as np

class Plot:
    def __init__(self, geometry, dem, physics):
        self.dts = []
        self.dx = geometry.dx
        self.elements = geometry.elements

        self.x_centres = geometry.centres()

        self.x_faces = geometry.interfaces()
        self.x_faces = np.repeat(self.x_faces, 2)
        self.x_faces = self.x_faces[1:-1]

        self.zs = dem.piecewise()
        self.velocity = physics.velocity

        plt.ion()
        plt.show()

    def __call__(self, state):
        Us = state.piecewise()

        plt.clf()
        plt.subplot(511)
        plt.gca().set_title('Elevation')
        plt.plot(self.x_faces, [z for z in self.zs], color='C2')
        plt.plot(self.x_faces, [z+U.h for z, U in zip(self.zs, Us)], color='C0')

        plt.subplot(512)
        plt.gca().set_title('Water depth')
        plt.plot(self.x_faces, [U.h for U in Us])
        plt.scatter(self.x_centres, [U.h.const for U in state], s=6)

        plt.subplot(513)
        plt.gca().set_title('Discharge')
        plt.plot(self.x_faces, [U.q for U in Us])
        plt.scatter(self.x_centres, [U.q.const for U in state], s=6)

        plt.subplot(514)
        plt.gca().set_title('Speed')
        plt.plot(self.x_centres, [abs(self.velocity(U.const())) for U in state])

        plt.subplot(515)
        plt.gca().set_title('dt')
        plt.yscale("log")
        plt.plot(*zip(*self.dts))

        plt.pause(1e-4)

    def block(self):
        plt.show(block=True)
