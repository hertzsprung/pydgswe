import matplotlib.pyplot as plt

def plot(xx, x_interface, z0, z1, h0, h1, q0, q1, tolh):
    plt.clf()
    plt.subplot(411)
    plt.gca().set_title("Elevation")
    plt.plot(xx, z0, color='C2')
    plt.plot(xx, [z+h for z, h in zip(z0, h0)], color='C0')
    plt.subplot(412)
    plt.gca().set_title("Depth")
    plt.plot(xx, h0)
    plt.subplot(413)
    plt.gca().set_title("Discharge")
    plt.plot(xx, q0)

    plt.subplot(414)
    plt.gca().set_title("Velocity")
    def velocity(h, q):
        return 0.0 if h <= tolh else q/h
    plt.plot(xx, [velocity(h, q) for h, q in zip(h0, q0)])

    plt.pause(1e-4)
